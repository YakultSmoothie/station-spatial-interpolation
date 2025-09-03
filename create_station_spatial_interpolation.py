#!/usr/bin/env python3
# ============================================================================================
# ==== INFORMATION ========
# ========================
# 檔名: create_station_spatial_interpolation.py
# 功能: 測站資料空間內插分析，支援多種內插方法與多時間點批量處理
# 建立日期: create 2025-09-03
# 更新日期: update 2025-09-03 (加入UTC時間轉換功能)
#
# Description:
#   此程式執行下列功能：
#   - 讀取測站基本資訊與觀測資料
#   - 執行空間內插分析(linear, cubic, nearest, idw, kriging)
#   - 支援當地時間轉UTC時間功能
#   - 輸出CSV、NetCDF格式結果檔案
#   - 產生科學期刊標準的視覺化圖表
#   - 支援多時間點批量處理
# ============================================================================================
import os
import time
from datetime import datetime, timedelta
start_time = time.time()
start_datetime = datetime.now()
print(f"{'='*80}")
print(f"Station Data Spatial Interpolation Analysis")
print(f"Started at: {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'='*80}")
#===========================================================================================
import pandas as pd
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')
import argparse
#===========================================================================================

def parse_arguments():
    """解析命令列參數"""
    parser = argparse.ArgumentParser(
        description='測站資料空間內插分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用範例:
  # 基本使用
  python3 create_station_spatial_interpolation.py -i ./data_station.txt -i2 ./data_pp01.txt

  # 變數, 時間, 目標網格 and 輸出路徑
  python3 create_station_spatial_interpolation.py -i ./data_station.txt -i2 ./data_pp01.txt -V PP01 -T 2025073001,2025073002,2025073003 -LL 120.0,122.12,21.9,25.4 -n ./output

  # 當地時間轉UTC時間 (台灣時間 UTC+8)
  python3 create_station_spatial_interpolation.py -i ./data_station.txt -i2 ./data_pp01.txt -V PP01 -T 2025073108 -utc +8

作者: CYC 2025-09-03 [v1.0]
        """)

    parser.add_argument('-i', '--input', required=True, help='測站基本資訊檔案路徑 (必要參數)')
    parser.add_argument('-i2', '--input2', required=True, help='觀測資料檔案路徑 (必要參數)')
    parser.add_argument('-V', '--variable', default='PP01', help='目標變數名稱 (預設: PP01)')
    parser.add_argument('-T', '--time', default='2025073001,2025073002', help='目標時間，逗號分隔 (預設: 2025073001,2025073002)')
    parser.add_argument('-LL', '--lonlat', default='120.0,122.12,21.9,25.4', help='網格經緯度範圍 lon_min,lon_max,lat_min,lat_max (預設: 120.0,122.12,21.9,25.4)')
    parser.add_argument('-n', '--output_dir', default='./output', help='輸出目錄路徑 (預設: ./output)')
    parser.add_argument('-r', '--resolution', type=float, default=0.03, help='網格解析度 (度) (預設: 0.03)')
    parser.add_argument('-m', '--method', choices=['linear', 'cubic', 'nearest', 'idw', 'kriging'], default='linear', help='內插方法 (預設: linear)')
    parser.add_argument('--idw_power', type=float, default=3.0, help='IDW方法的冪次指數 (預設: 3.0)')
    parser.add_argument('--max_distance', type=float, default=50.0, help='IDW方法的最大搜索距離 (公里) (預設: 50.0)')
    parser.add_argument('-nd', '--no_vis', action='store_true', help='不產生視覺化圖表')
    parser.add_argument('-utc', '--utc_offset', type=str, default=None, 
                       help='當地時間轉UTC時間偏移量，格式: +N 或 -N (小時)，如: +8 表示UTC+8轉UTC')

    return parser.parse_args()

# -----------------
# PARAMETERS - 可控變數設定
# -----------------
def setup_parameters(args):
    """根據輸入參數設定全域變數"""
    global IFN1, IFN2, TARGET_VARIABLE, TARGET_TIME_ORI, LON_RANGE, LAT_RANGE
    global GRID_RESOLUTION, INTERPOLATION_METHOD, IDW_POWER, MAX_DISTANCE, OUTPUT_DIR, RUN_CIV
    global UTC_OFFSET, CONVERT_TO_UTC
    
    # 必要輸入參數
    IFN1 = args.input                    # 測站基本資訊檔案  required arguments
    IFN2 = args.input2                   # 觀測資料檔案  required arguments
    TARGET_VARIABLE = args.variable      # 目標變數名稱  required arguments
    TARGET_TIME_ORI = args.time          # 目標時間  required arguments
    
    # UTC轉換參數處理
    if args.utc_offset:
        CONVERT_TO_UTC = True
        try:
            UTC_OFFSET = int(args.utc_offset)
            print(f"UTC conversion enabled: Local time UTC{args.utc_offset} -> UTC")
        except ValueError:
            raise ValueError(f"UTC offset format should be +N or -N (hours), got: {args.utc_offset}")
    else:
        CONVERT_TO_UTC = False
        UTC_OFFSET = 0
    
    # 解析經緯度範圍
    lonlat_parts = [float(x) for x in args.lonlat.split(',')]
    if len(lonlat_parts) != 4:
        raise ValueError("經緯度範圍格式應為: lon_min,lon_max,lat_min,lat_max")
    LON_RANGE = [lonlat_parts[0], lonlat_parts[1]]    # 網格經度範圍 [最小值, 最大值]
    LAT_RANGE = [lonlat_parts[2], lonlat_parts[3]]    # 網格緯度範圍 [最小值, 最大值]

    # 選用參數
    GRID_RESOLUTION = args.resolution       # 網格間距 (度)
    INTERPOLATION_METHOD = args.method      # 內插方法: 'linear', 'cubic', 'nearest', 'idw', 'kriging'
    IDW_POWER = args.idw_power             # 距離的冪次指數 for idw only. 
    MAX_DISTANCE = args.max_distance       # 最大搜索距離 (km)  for idw only.
    OUTPUT_DIR = args.output_dir
    # 是否視覺化
    RUN_CIV = not args.no_vis              # key setting create_interpolation_visualization

    # 確保輸出目錄存在
    os.makedirs(OUTPUT_DIR, exist_ok=True)

# -----------------
# UTC TIME CONVERSION - UTC時間轉換功能
# -----------------
def convert_local_to_utc_time(local_time_str, utc_offset):
    """將當地時間轉換為UTC時間
    
    Args:
        local_time_str (str): 當地時間格式 (YYYYMMDDHH)
        utc_offset (int): UTC偏移量，正數表示東時區，負數表示西時區
    
    Returns:
        str: UTC時間格式 (YYYYMMDDHH)
    
    Examples:
        convert_local_to_utc_time('2025073108', 8) -> '2025073100'  # UTC+8 -> UTC
        convert_local_to_utc_time('2025073108', -5) -> '2025073113' # UTC-5 -> UTC
    """
    print(f"Converting local time to UTC:")
    print(f"    Input local time: {local_time_str} (UTC{utc_offset:+d})")
    
    # 解析當地時間
    local_datetime = pd.to_datetime(local_time_str, format='%Y%m%d%H')
    
    # 轉換為UTC時間（減去時區偏移）
    utc_datetime = local_datetime - timedelta(hours=utc_offset)
    
    # 轉換為十碼格式
    utc_time_str = utc_datetime.strftime('%Y%m%d%H')
    
    print(f"    Output UTC time: {utc_time_str}")
    print(f"    Conversion: {local_datetime.strftime('%Y-%m-%d %H:%M')} (Local) -> {utc_datetime.strftime('%Y-%m-%d %H:%M')} (UTC)")
    
    return utc_time_str

def process_time_with_utc_conversion(time_ori, ctu):
    """處理時間轉換，包含UTC轉換功能
    
    Args:
        time_ori (str): 原始時間格式 (YYYYMMDDHH)
    
    Returns:
        tuple: (處理後的時間字串, 用於檔案命名的時間字串)
    """
    if ctu:
        # 當地時間轉UTC
        utc_time_str = convert_local_to_utc_time(time_ori, UTC_OFFSET)
        # 轉換為標準時間格式
        target_time = change_to_new_time(time_ori)
        # 用於檔案命名的時間字串（包含UTC標記）
        filename_time = f"{utc_time_str}_UTC"
        print(f"    Final target time for analysis: {target_time}")
        return target_time, filename_time
    else:
        # 不進行UTC轉換
        target_time = change_to_new_time(time_ori)
        filename_time = time_ori
        return target_time, filename_time

# -----------------
# OPEN - 讀取測站資料
# -----------------
def change_to_new_time(time_ori):
    """將十碼時間格式轉換為標準時間格式
    
    Args:
        time_ori (str): 十碼時間格式 (YYYYMMDDHH)
    
    Returns:
        str: 標準時間格式 (YYYY-MM-DD HH:00:00)
    """
    print(f"Converting time format from {time_ori}")
    
    # 確保輸入是字串格式
    time_str = str(time_ori)
    
    # 檢查格式是否正確
    if len(time_str) != 10:
        raise ValueError(f"Time format should be 10 digits (YYYYMMDDHH), got: {time_str}")
    
    # 解析各個時間組件
    year = time_str[:4]
    month = time_str[4:6]
    day = time_str[6:8]
    hour = time_str[8:10]
    
    # 處理24小時的情況
    if hour == "24":
        # 將24小時轉換為次日00小時
        from datetime import datetime, timedelta
        base_date = datetime(int(year), int(month), int(day))
        next_day = base_date + timedelta(days=1)
        converted_time = next_day.strftime("%Y-%m-%d 00:00:00")
        print(f"    24-hour converted to next day: {converted_time}")
    else:
        converted_time = f"{year}-{month}-{day} {hour}:00:00"
        print(f"    Converted to: {converted_time}")
    
    return converted_time

def change_to_old_time(time_standard):
    """將標準時間格式轉換為十碼時間格式
    
    Args:
        time_standard (str): 標準時間格式 (YYYY-MM-DD HH:00:00)
    
    Returns:
        str: 十碼時間格式 (YYYYMMDDHH)
    """
    print(f"Converting time format from {time_standard}")
    
    # 確保輸入是字串格式
    time_str = str(time_standard).strip()
    
    # 處理可能的時間格式變化
    if 'T' in time_str:
        # ISO格式: YYYY-MM-DDTHH:MM:SS
        date_part, time_part = time_str.split('T')
        hour = time_part.split(':')[0]
    elif ' ' in time_str:
        # 標準格式: YYYY-MM-DD HH:MM:SS
        date_part, time_part = time_str.split(' ')
        hour = time_part.split(':')[0]
    else:
        # 只有日期: YYYY-MM-DD
        date_part = time_str
        hour = "00"
    
    # 解析日期部分
    try:
        year, month, day = date_part.split('-')
    except ValueError:
        raise ValueError(f"Date format should be YYYY-MM-DD, got: {date_part}")
    
    # 驗證各個時間組件
    if len(year) != 4 or not year.isdigit():
        raise ValueError(f"Year should be 4 digits, got: {year}")
    if len(month) != 2 or not month.isdigit() or not (1 <= int(month) <= 12):
        raise ValueError(f"Month should be 2 digits (01-12), got: {month}")
    if len(day) != 2 or not day.isdigit() or not (1 <= int(day) <= 31):
        raise ValueError(f"Day should be 2 digits (01-31), got: {day}")
    if len(hour) != 2 or not hour.isdigit() or not (0 <= int(hour) <= 23):
        raise ValueError(f"Hour should be 2 digits (00-23), got: {hour}")
    
    converted_time = f"{year}{month}{day}{hour}"
    print(f"    Converted to: {converted_time}")
    
    return converted_time

def load_station_info(station_file):
    """讀取測站基本資訊檔案"""
    print(f"Step 1A: Reading station information from {station_file}")
    
    try:
        # 嘗試不同的分隔符號
        df_station = None
        separator_used = None
        
        # 先嘗試空白分隔符
        try:
            df_station = pd.read_csv(station_file, sep=r'\s+')
            separator_used = "whitespace"
            print(f"    Using whitespace separator")
        except:
            try:
                df_station = pd.read_csv(station_file, sep='|')
                separator_used = "pipe"
                print(f"    Using pipe separator (|)")
            except:
                try:
                    df_station = pd.read_csv(station_file, sep=',')
                    separator_used = "comma"
                    print(f"    Using comma separator")
                except:
                    raise Exception("Unable to parse file with any common separator")
        
        print(f"    Station file loaded successfully:")
        print(f"    Columns: {df_station.columns.tolist()}")
        print(f"    Number of stations: {len(df_station)}")
        
        # 檢查欄位數量
        if len(df_station.columns) != 4:
            print(f"    Warning: Expected 4 columns, but found {len(df_station.columns)}")
            if len(df_station.columns) < 4:
                print(f"    Error: Insufficient columns for station data")
                return None
            
        # 檢查是否已有適當的欄位名稱，如果沒有才重新命名
        expected_cols = ['station_id', 'longitude', 'latitude', 'elevation']
        if list(df_station.columns) != expected_cols:
            print(f"    Renaming columns to standard format")
            df_station.columns = expected_cols[:len(df_station.columns)]
            
        print(f"    Station coordinate range:")
        print(f"        Longitude: {df_station['longitude'].min():.3f} to {df_station['longitude'].max():.3f}")
        print(f"        Latitude: {df_station['latitude'].min():.3f} to {df_station['latitude'].max():.3f}")
        
        return df_station
        
    except Exception as e:
        print(f"    Error reading station file: {e}")
        return None

def load_observation_data(data_file):
    """
    讀取觀測資料檔案
    若觀測資料有紀錄小時為24的自動改為後一天的00
    """
    print(f" ")
    print(f"Step 1B: Reading observation data from {data_file}")
    
    try:
        # 嘗試不同的分隔符號
        try:
            df_obs = pd.read_csv(data_file, sep=r'\s+')
            print(f"    Using whitespace separator")
        except:
            try:
                df_obs = pd.read_csv(data_file, sep='|')
                print(f"    Using pipe separator (|)")
            except:
                df_obs = pd.read_csv(data_file, sep=',')
                print(f"    Using comma separator")
        
        print(f"    Observation file loaded successfully:")
        print(f"    Columns: {df_obs.columns.tolist()}")
        print(f"    Number of records: {len(df_obs)}")
        
        # 檢查是否有24小時的資料
        hours = df_obs['yyyymmddhh'] % 100
        count_24h = (hours > 23).sum()
        if count_24h > 0:
            print(f"    Found {count_24h} records with 24-hour format (will be converted to next day 00:00)")
        
        print(f"    Converting time format...")
        
        # 向量化處理 - 更快的方式
        # 先將所有24小時改為00小時
        df_temp = df_obs['yyyymmddhh'].copy()
        mask_24h = (hours == 24)
        df_temp[mask_24h] = df_temp[mask_24h] - 24
        
        # 轉換為datetime
        df_obs['datetime'] = pd.to_datetime(df_temp.astype(str), format='%Y%m%d%H')
        
        # 對24小時的資料加1天
        df_obs.loc[mask_24h, 'datetime'] += pd.Timedelta(days=1)
        
        # 保留原始欄位並重新命名
        df_obs = df_obs.rename(columns={'yyyymmddhh': 'original_yyyymmddhh'})
        df_obs['yyyymmddhh'] = df_obs['datetime']
        
        print(f"    Time range: {df_obs['yyyymmddhh'].min()} to {df_obs['yyyymmddhh'].max()}")
        print(f"    Available variables: {[col for col in df_obs.columns if col not in ['stno', 'yyyymmddhh', 'original_yyyymmddhh', 'datetime', 'station_id']]}")
        
        # 標準化站號欄位名稱
        if 'stno' in df_obs.columns:
            df_obs = df_obs.rename(columns={'stno': 'station_id'})
            print(f"    Renamed 'stno' to 'station_id' for consistency")
        
        return df_obs
        
    except Exception as e:
        print(f"    Error reading observation file: {e}")
        return None

# -----------------
# DEFINE - 資料整合與前處理
# -----------------
def integrate_station_data(df_station, df_obs, target_variable, target_time):
    """Step 1: 資料整合與前處理"""
    print(f"\nStep 2: Integrating station and observation data")
    print(f"    Target variable: {target_variable}")
    print(f"    Target time: {target_time}")
    
    # 轉換目標時間為datetime格式
    target_datetime = pd.to_datetime(target_time)
    
    # 篩選指定時間的資料
    df_obs_filtered = df_obs[df_obs['yyyymmddhh'] == target_datetime].copy()
    print(f"    Records at target time: {len(df_obs_filtered)}")
    
    # 檢查是否有時間上的資料可用
    if len(df_obs_filtered) == 0:
        print(f"    Warning: No data found at target time {target_time}")
        print(f"        Available times: {sorted(df_obs['yyyymmddhh'].unique())}")
        return None
    
    # 檢查測站匹配情況
    # breakpoint()
    obs_stations = set(df_obs_filtered['station_id'].unique())
    station_stations = set(df_station['station_id'].unique())
    missing_stations = obs_stations - station_stations
    if missing_stations:
        print(f"    Warning: {len(missing_stations)} stations in observation data not found in station info:")
        print(f"        Missing stations: {sorted(list(missing_stations))}")
    
    # 檢查目標變數是否存在
    if target_variable not in df_obs.columns:
        print(f"    Error: Variable '{target_variable}' not found in observation data")
        print(f"        Available variables: {df_obs.columns[2:].tolist()}")
        return None
    
    # 標準化觀測資料欄位名稱
    obs_cols = ['station_id', 'datetime'] + df_obs.columns[2:].tolist()
    df_obs_filtered.columns = obs_cols
    
    # 合併測站資訊與觀測資料 只有同時存在於兩個資料框中的測站才會被保留在最終結果中。
    df_integrated = pd.merge(df_station, df_obs_filtered, on='station_id', how='inner')  # 進行內部合併
    
    # 移除含有缺失值的記錄
    initial_count = len(df_integrated)
    df_integrated = df_integrated.dropna(subset=[target_variable, 'longitude', 'latitude'])
    final_count = len(df_integrated)
    
    print(f"    Integration completed:")
    print(f"        Initial records: {initial_count}")
    print(f"        Valid records after removing NaN: {final_count}")
    print(f"        Variable '{target_variable}' range: {df_integrated[target_variable].min():.2f} to {df_integrated[target_variable].max():.2f}")
    
    # 篩選在目標區域內的測站
    # lon_mask = (df_integrated['longitude'] >= LON_RANGE[0]) & (df_integrated['longitude'] <= LON_RANGE[1])
    # lat_mask = (df_integrated['latitude'] >= LAT_RANGE[0]) & (df_integrated['latitude'] <= LAT_RANGE[1])
    # df_integrated = df_integrated[lon_mask & lat_mask]

    # 只保留必要欄位
    keep_cols = ['station_id', 'longitude', 'latitude', 'elevation', 'datetime', target_variable]
    df_integrated = df_integrated[keep_cols]
    
    print(f"        Records within target domain: {len(df_integrated)}")
    
    return df_integrated

# -----------------
# DEFINE - 建立目標網格系統
# -----------------
def create_target_grid(lon_range, lat_range, resolution):
    """Step 2: 建立目標網格系統"""
    print(f"\nStep 1C: Creating target grid system")
    print(f"    Longitude range: {lon_range}")
    print(f"    Latitude range: {lat_range}")
    print(f"    Grid resolution: {resolution} degrees")
    
    # 建立一維網格座標
    lon_grid = np.arange(lon_range[0], lon_range[1] + resolution, resolution)
    lat_grid = np.arange(lat_range[0], lat_range[1] + resolution, resolution)
    
    # 建立二維網格座標
    lon_mesh, lat_mesh = np.meshgrid(lon_grid, lat_grid)
    
    print(f"    Grid dimensions: {len(lat_grid)} x {len(lon_grid)} (lat x lon)")
    print(f"    Total grid points: {lon_mesh.size}")
    
    grid_info = {
        'lon_grid': lon_grid,
        'lat_grid': lat_grid,
        'lon_mesh': lon_mesh,
        'lat_mesh': lat_mesh
    }
    
    return grid_info

# -----------------
# DEFINE - 執行空間內插
# -----------------
def perform_spatial_interpolation(df_integrated, grid_info, target_variable, method='linear'):
    """Step 3: 執行空間內插"""
    print(f"\nStep 3: Performing spatial interpolation")
    print(f"    Interpolation method: {method}")
    print(f"    Number of observation points: {len(df_integrated)}")
    
    # 提取測站座標與觀測值
    station_coords = df_integrated[['longitude', 'latitude']].values
    station_values = df_integrated[target_variable].values
    
    # 提取網格座標
    grid_coords = np.column_stack([grid_info['lon_mesh'].ravel(), 
                                   grid_info['lat_mesh'].ravel()])
    
    # 執行內插
    try:
        if method == 'idw':
            # 使用反距離權重法 (Inverse Distance Weighting)
            from scipy.spatial.distance import cdist
            print(f"    Using Inverse Distance Weighting (IDW) method")
            interpolated_values = idw_interpolation(station_coords, station_values, 
                                                   grid_coords, power=IDW_POWER, max_distance=MAX_DISTANCE)
            #do more time, and use max
            interpolated_grid = interpolated_values.reshape(grid_info['lat_mesh'].shape)  # 重塑為二維網格
            
        elif method == 'kriging':
            # 使用克利金內插法
            from pykrige.ok import OrdinaryKriging
            print(f"    Using Ordinary Kriging method")
                
            # 經緯度座標
            lons = station_coords[:, 0]
            lats = station_coords[:, 1]
            values = station_values
                         
            # 建立最佳Kriging模型
            OK = OrdinaryKriging(lons, lats, values, variogram_model='spherical', enable_plotting=False, verbose=False, coordinates_type='geographic')
            
            # 在網格上內插
            grid_lon = grid_info['lon_mesh'][0, :]
            grid_lat = grid_info['lat_mesh'][:, 0]
            interpolated_grid, uncertainty = OK.execute('grid', grid_lon, grid_lat)
            
            print(f"        Interpolated value range: {np.nanmin(interpolated_grid):.2f} to {np.nanmax(interpolated_grid):.2f}")
            print(f"        Kriging variance range: {np.nanmin(uncertainty):.4f} to {np.nanmax(uncertainty):.4f}")

        else:
            # 使用scipy的griddata方法
            from scipy.interpolate import griddata
            print(f"    Using scipy griddata with method: {method}")
            interpolated_values = griddata(station_coords, station_values, 
                                         grid_coords, method=method, 
                                         fill_value=np.nan)
            interpolated_grid = interpolated_values.reshape(grid_info['lat_mesh'].shape)  # 重塑為二維網格            
        
        # 統計內插結果
        valid_points = ~np.isnan(interpolated_grid)
        valid_count = np.sum(valid_points)
        total_count = interpolated_grid.size
        
        print(f"    Interpolation completed:")
        print(f"        Valid interpolated points: {valid_count} / {total_count} ({valid_count/total_count*100:.1f}%)")
        print(f"        Interpolated value range: {np.nanmin(interpolated_grid):.2f} to {np.nanmax(interpolated_grid):.2f}")
        
        return interpolated_grid
        
    except Exception as e:
        print(f"    Error during interpolation: {e}")
        return None

def calculate_geographic_distances_fully_vectorized(target_coords, known_coords):
    """
    完全向量化的地理距離計算
    使用 Haversine 公式
    """
    # 轉換為弧度
    target_lat_rad = np.radians(target_coords[:, 1:2])  # shape: (n_target, 1)
    target_lon_rad = np.radians(target_coords[:, 0:1])  # shape: (n_target, 1)
    known_lat_rad = np.radians(known_coords[np.newaxis, :, 1])  # shape: (1, n_known)
    known_lon_rad = np.radians(known_coords[np.newaxis, :, 0])  # shape: (1, n_known)
    
    # Haversine 公式
    dlat = known_lat_rad - target_lat_rad  # broadcasting to (n_target, n_known)
    dlon = known_lon_rad - target_lon_rad  # broadcasting to (n_target, n_known)
    
    a = (np.sin(dlat/2)**2 + 
         np.cos(target_lat_rad) * np.cos(known_lat_rad) * np.sin(dlon/2)**2)
    
    c = 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))  # clip to avoid numerical errors
    
    # 地球半徑（公里）
    R = 6371.009
    distances_km = R * c
    
    return distances_km

def idw_interpolation(known_coords, known_values, target_coords, power=2, max_distance=None):
    """反距離權重法內插"""
    # 計算距離矩陣

    # 座標距離(歐幾里得距離)
    # distances_deg = cdist(target_coords, known_coords) 
    # 輸出距離資訊
    # print(f"    Distance matrix information (deg):")
    # print(f"        Shape: {distances_deg.shape} (target points × observation points)")
    # print(f"        Distance range: {np.min(distances_deg):.6f} to {np.max(distances_deg):.6f} degrees")
    # print(f"        Mean distance: {np.mean(distances_deg):.6f} degrees")
# 
    # 計算地理距離（公里）
    distances_km = calculate_geographic_distances_fully_vectorized(target_coords, known_coords)

    # 使用地理距離進行IDW計算
    distances = distances_km.copy()

    # 輸出距離資訊
    print(f"    Distance matrix information (km)")
    print(f"        Shape: {distances.shape} (target points × observation points)")
    print(f"        Distance range: {np.min(distances):.6f} to {np.max(distances):.6f} degrees")
    print(f"        Mean distance: {np.mean(distances):.6f} degrees")
    
    # 應用最大距離限制
    if max_distance is not None:
        distances[distances > max_distance] = np.inf
        print(f"        Applying max distance limit: {max_distance:.1f} km")
        excluded_points = np.sum(distances == np.inf)
        #print(f"        Excluded {excluded_points} distance pairs due to max distance limit")
    
    # 避免除零錯誤
    distances[distances == 0] = 1e-10
    
    # 計算權重 (距離的負幾次方)
    print(f"        Computing IDW weights with power parameter: {power}")
    weights = 1.0 / (distances ** power)
    print(f"        Weight range: {np.min(weights[weights > 0]):.6f} to {np.max(weights[~np.isinf(weights)]):.6f}")

    # 正規化權重
    weights_sum = np.sum(weights, axis=1, keepdims=True)
    weights_normalized = weights / weights_sum
    
    # 計算加權平均
    interpolated_values = np.sum(weights_normalized * known_values[np.newaxis, :], axis=1)
    
    return interpolated_values

# -----------------
# WRITE AND SAVE - 輸出結果檔案
# -----------------
def save_integrated_dataframe(df_integrated, output_dir, target_variable, filename_time):
    """保存整合資料框為CSV檔案"""
    print(f"\nStep 4A: Saving integrated dataframe")
    
    output_dir_csv = f'{output_dir}/csv'
    os.makedirs(output_dir_csv, exist_ok=True)
    output_file = os.path.join(output_dir_csv, f'integrated_station_data_{target_variable}_{filename_time}.csv')
    df_integrated.to_csv(output_file, index=False)
    
    print(f"    Integrated dataframe saved to: {output_file}")
    print(f"    File size: {os.path.getsize(output_file)/1024:.1f} KB")
    
    return output_file

def save_gridded_data(interpolated_grid, grid_info, target_variable, target_time, filename_time, output_dir, CONVERT_TO_UTC):
    """保存網格資料為NetCDF檔案"""
    print(f"\nStep 4B: Saving gridded data (to .nc)")
   
    # 建立xarray Dataset
    ds = xr.Dataset(
        data_vars={
            target_variable: (['lat', 'lon'], interpolated_grid)
        },
        coords={
            'lon': grid_info['lon_grid'],
            'lat': grid_info['lat_grid'],
            'time': pd.to_datetime(target_time)
        },
        attrs={
            'title': f'Interpolated {target_variable} data',
            'description': f'Spatial interpolation of station observations',
            'target_time': target_time,
            'interpolation_method': INTERPOLATION_METHOD,
            'grid_resolution': GRID_RESOLUTION,
            'creation_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'Conventions': 'CF-1.6',
            'history': f'Created by spatial interpolation on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
            'utc_conversion': f'UTC conversion applied: {CONVERT_TO_UTC}' if CONVERT_TO_UTC else 'No UTC conversion applied',
            'utc_offset': f'UTC{UTC_OFFSET:+d}' if CONVERT_TO_UTC else 'N/A'
        }
    )
    
    # 設定座標屬性
    ds.lon.attrs = {
        'long_name': 'longitude',
        'standard_name': 'longitude',
        'units': 'degrees_east',
        'axis': 'X'
    }
    
    ds.lat.attrs = {
        'long_name': 'latitude', 
        'standard_name': 'latitude',
        'units': 'degrees_north',
        'axis': 'Y'
    }
    
    ds.time.attrs = {
        'long_name': 'time',
        'standard_name': 'time',
        'axis': 'T',
    }
    
    ds[target_variable].attrs = {
        'long_name': f'Interpolated {target_variable}',
        'standard_name': target_variable.lower().replace(' ', '_'),
        'units': 'unknown',
        'interpolation_method': INTERPOLATION_METHOD,
        'coordinates': 'lat lon time'
    }
    
    # 分別轉換座標和資料變數的資料類型
    ds = ds.assign_coords({
        'lat': ds.lat.astype('float32'),
        'lon': ds.lon.astype('float32')
    })
    ds[target_variable] = ds[target_variable].astype('float32')
    
    # 定義 encoding 參數 (不壓縮)
    encoding = {
        'lat': {'dtype': 'float32', '_FillValue': None},
        'lon': {'dtype': 'float32', '_FillValue': None},
        'time': {'dtype': 'float64', '_FillValue': None},
        target_variable: {
            'dtype': 'float32',
            '_FillValue': np.nan
        }
    }
    
    # 保存為NetCDF
    output_dir_nc = f'{output_dir}/nc'
    os.makedirs(output_dir_nc, exist_ok=True)
    output_file = os.path.join(output_dir_nc, f'interpolated_{target_variable}_{filename_time}.nc')
    ds.to_netcdf(output_file, format='NETCDF4', encoding=encoding)
    
    print(f"    Gridded data saved to: {output_file}")
    print(f"    File size: {os.path.getsize(output_file)/1024:.1f} KB")
    print(f"    Dataset coordinates:")
    print(f"        Longitude: {ds.lon.min().values:.3f} to {ds.lon.max().values:.3f} degrees_east")
    print(f"        Latitude: {ds.lat.min().values:.3f} to {ds.lat.max().values:.3f} degrees_north")
    print(f"        Time: {ds.time.values} ")
    print(f"    GrADS compatibility: ENABLED")
    
    return output_file, ds

# -----------------
# PLOT AND SAVE - 視覺化結果
# -----------------
def create_interpolation_visualization(df_integrated, grid_info, interpolated_grid, 
                                      target_variable, target_time, filename_time, output_dir, CONVERT_TO_UTC):
    """創建內插結果視覺化圖表"""
    print(f"\nCreating interpolation visualization")
    import matplotlib.pyplot as plt
    import matplotlib.colors as clr
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    # 設置字型大小和圖形樣式
    FONT_SIZE = 14
    plt.rcParams.update({
       'font.size': FONT_SIZE,
       'axes.titlesize': FONT_SIZE * 1.2,
       'axes.labelsize': FONT_SIZE * 1.0,  
       'xtick.labelsize': FONT_SIZE * 0.9,
       'ytick.labelsize': FONT_SIZE * 0.9,
       'legend.fontsize': FONT_SIZE * 0.9,
    })

    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), 
                                   subplot_kw={'projection': ccrs.PlateCarree()})
    
    # 設定地圖範圍
    extent = [LON_RANGE[0], LON_RANGE[1], LAT_RANGE[0], LAT_RANGE[1]]
    
    for ax in [ax1, ax2]:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, linewidth=1)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.LAND, alpha=0.3)
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5, zorder=2)
        gl.xlabel_style = {'size': FONT_SIZE * 0.6}  # 調整經度標籤字體大小
        gl.ylabel_style = {'size': FONT_SIZE * 0.6}  # 調整緯度標籤字體大小
        gl.top_labels = False      # 不在上方顯示經度標籤
        gl.right_labels = False    # 不在右邊顯示緯度標籤

    # set color：
    vmin = df_integrated[target_variable].min()
    vmax = df_integrated[target_variable].max()
    # vmin2 = np.min(interpolated_grid)
    # vmax2 = np.max(interpolated_grid)
    # vmin = min(vmin, vmin2)
    # vmax = max(vmax, vmax2)


    # 建立一個共用的 norm for color
    norm = clr.Normalize(vmin=vmin, vmax=vmax)
    
    # 左圖：scatter 使用 norm
    scatter = ax1.scatter(df_integrated['longitude'], df_integrated['latitude'], 
                         c=df_integrated[target_variable], s=25, 
                         cmap='RdYlBu_r', edgecolors='black', linewidth=0.3,
                         transform=ccrs.PlateCarree(), zorder=5,
                         norm=norm)
    
    # 根據是否有UTC轉換來設定標題
    time_label = f"{target_time}"
    ax1.set_title(f'Station Observations\n{target_variable} at {time_label}', 
                 fontsize=FONT_SIZE * 1.1)
    
    # 右圖：contourf 使用 norm
    contour = ax2.pcolormesh(grid_info['lon_mesh'], grid_info['lat_mesh'], 
                          interpolated_grid, cmap='RdYlBu_r', 
                          transform=ccrs.PlateCarree(), zorder=1,
                          norm=norm)
    
    # 疊加測站位置
    ax2.scatter(df_integrated['longitude'], df_integrated['latitude'], 
               c='black', s=3, marker='o', 
               transform=ccrs.PlateCarree(), zorder=5)
    
    ax2.set_title(f'Interpolated Grid\nMethod: {INTERPOLATION_METHOD}', 
                  fontsize=FONT_SIZE * 1.1)
    
    # 添加色標各自的 colorbar
    cbar1 = plt.colorbar(scatter, ax=ax1, orientation='horizontal', extend='both',
                        pad=0.05, aspect=40, shrink=0.8)
    cbar1.set_label(f'{target_variable}', fontsize=FONT_SIZE)

    cbar2 = plt.colorbar(contour, ax=ax2, orientation='horizontal', extend='both',
                        pad=0.05, aspect=40, shrink=0.8)
    cbar2.set_label(f'{target_variable}', fontsize=FONT_SIZE)

    
    # 加粗外框
    for ax in [ax1, ax2]:
        for spine in ax.spines.values():
            spine.set_linewidth(3)
            spine.set_zorder(99)
    
    plt.tight_layout()
    
    # 保存圖片
    output_dir_png = f'{output_dir}/png'
    os.makedirs(output_dir_png, exist_ok=True)
    output_file = os.path.join(output_dir_png, f'interpolation_visualization_{target_variable}_{filename_time}_{INTERPOLATION_METHOD}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"    Visualization saved to: {output_file}")
    
    return output_file

# -----------------
# MAIN - 主執行流程
# -----------------
def main():
    """主執行流程 - 支援多時間點處理與UTC轉換"""
    
    # 解析命令列參數
    args = parse_arguments()
    
    # 設定全域參數
    setup_parameters(args)
    
    # 處理多時間點輸入
    if isinstance(TARGET_TIME_ORI, str) and ',' in TARGET_TIME_ORI:
        TARGET_TIME_ORIs = [time.strip() for time in TARGET_TIME_ORI.split(',')]
    elif isinstance(TARGET_TIME_ORI, str):
        TARGET_TIME_ORIs = [TARGET_TIME_ORI]
    else:
        TARGET_TIME_ORIs = TARGET_TIME_ORI if isinstance(TARGET_TIME_ORI, list) else [TARGET_TIME_ORI]
    
    print(f"{'='*80}")
    print(f"Processing parameters:")
    print(f"    Station file: {IFN1}")
    print(f"    Data file: {IFN2}")
    print(f"    Target variable: {TARGET_VARIABLE}")
    print(f"    Target times: {TARGET_TIME_ORIs}")
    print(f"    UTC conversion: {'Yes' if CONVERT_TO_UTC else 'No'}")
    if CONVERT_TO_UTC:
        print(f"    UTC offset: UTC{UTC_OFFSET:+d} -> UTC")
    print(f"    Longitude range: {LON_RANGE}")
    print(f"    Latitude range: {LAT_RANGE}")
    print(f"    Grid resolution: {GRID_RESOLUTION}°")
    print(f"    Interpolation method: {INTERPOLATION_METHOD}")
    print(f"    Output directory: {OUTPUT_DIR}")
    print(f"{'='*80}")
    
    # Step 0: 載入基礎資料 (只需載入一次)
    print(" ")
    df_station = load_station_info(IFN1)
    if df_station is None:
        return
    
    df_obs = load_observation_data(IFN2)
    if df_obs is None:
        return
    
    # 建立目標網格系統 (只需建立一次)
    print(" ")
    grid_info = create_target_grid(LON_RANGE, LAT_RANGE, GRID_RESOLUTION)
    print(f"    Grid system = target interpolation mesh:")
    print(f"        Longitude points: {len(grid_info['lon_grid'])}")
    print(f"        Latitude points: {len(grid_info['lat_grid'])}")
    
    # 儲存所有時間點的結果
    all_results = {}
    all_csv_files = []
    all_nc_files = []
    all_vis_files = []
    
    # 對每個時間點進行處理
    for i, TARGET_TIME_ORI_single in enumerate(TARGET_TIME_ORIs):
        print(f"\n{'='*60}")
        print(f"Processing time point {i+1}/{len(TARGET_TIME_ORIs)}: {TARGET_TIME_ORI_single}")
        print(f"{'='*60}")
        
        # 時間格式轉換（包含UTC轉換）
        TARGET_TIME, FILENAME_TIME = process_time_with_utc_conversion(TARGET_TIME_ORI_single, CONVERT_TO_UTC)
        print(f"    Target time for analysis: {TARGET_TIME}")
        print(f"    Filename identifier: {FILENAME_TIME}")
        
        # Step 1: 資料整合與前處理
        df_integrated = integrate_station_data(df_station, df_obs, TARGET_VARIABLE, TARGET_TIME)
        if df_integrated is None or len(df_integrated) == 0:
            print(f"Warning: No valid integrated data available for time {TARGET_TIME_ORI_single}")
            continue
        
        print(f"    Final integrated dataset = climatological station obs at target time:")
        print(f"        Columns: {df_integrated.columns.tolist()}")
        print(f"        Shape: {df_integrated.shape}")
        
        # Step 3: 執行空間內插
        print(" ")
        interpolated_grid = perform_spatial_interpolation(df_integrated, grid_info, 
                                                         TARGET_VARIABLE, INTERPOLATION_METHOD)
        if interpolated_grid is None:
            print(f"Error: Interpolation failed for time {TARGET_TIME_ORI_single}")
            continue
        
        print(f"    Interpolated grid = spatial analysis result:")
        print(f"        Grid shape: {interpolated_grid.shape}")
        print(f"        Value statistics: min={np.nanmin(interpolated_grid):.2f}, max={np.nanmax(interpolated_grid):.2f}, mean={np.nanmean(interpolated_grid):.2f}")
        
        # Step 4: 輸出結果檔案
        print(" ")
        csv_file = save_integrated_dataframe(df_integrated, OUTPUT_DIR, TARGET_VARIABLE, FILENAME_TIME)
        nc_file, ds = save_gridded_data(interpolated_grid, grid_info, TARGET_VARIABLE, 
                                       TARGET_TIME, FILENAME_TIME, OUTPUT_DIR, CONVERT_TO_UTC)
        
        # 創建視覺化
        vis_file = None
        if RUN_CIV:
            vis_file = create_interpolation_visualization(df_integrated, grid_info, interpolated_grid,
                                                         TARGET_VARIABLE, TARGET_TIME, FILENAME_TIME, OUTPUT_DIR, CONVERT_TO_UTC)
        
        # 儲存此時間點的結果
        all_results[TARGET_TIME_ORI_single] = {
            'df_integrated': df_integrated,
            'interpolated_grid': interpolated_grid,
            'dataset': ds,
            'target_time': TARGET_TIME,
            'filename_time': FILENAME_TIME
        }
        all_csv_files.append(csv_file)
        all_nc_files.append(nc_file)
        if vis_file:
            all_vis_files.append(vis_file)
        
        print(f"\nTime point {TARGET_TIME_ORI_single} completed successfully!")
        if CONVERT_TO_UTC:
            print(f"    Input local time: {TARGET_TIME_ORI_single} (UTC{UTC_OFFSET:+d})")
            print(f"    Converted to UTC: {FILENAME_TIME.replace('_UTC', '')}")
        print(f"    CSV: {csv_file}")
        print(f"    NetCDF: {nc_file}")
        if vis_file:
            print(f"    Visualization: {vis_file}")
    
    # 最終總結
    end_time = time.time()
    end_datetime = datetime.now()
    print(f"\n{'='*80}")
    print(f"Multi-time analysis completed successfully!")
    print(f"Processed {len(all_results)} time points: {list(all_results.keys())}")
    if CONVERT_TO_UTC:
        print(f"UTC conversion applied: Local time UTC{UTC_OFFSET:+d} -> UTC")
    print(f"Total output files created:")
    print(f"    CSV files: {len(all_csv_files)}")
    print(f"    NetCDF files: {len(all_nc_files)}")
    print(f"    Visualization files: {len(all_vis_files)}")
    print(f"{'-'*80}")
    print(f"End at: {end_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total processing time: {end_time - start_time:.2f} seconds")
    print(f"{'='*80}")
    
    return all_results, grid_info

if __name__ == "__main__":
    all_results, grid_info = main()
    
    # 如果需要檢查特定時間點的結果，可以這樣存取：
    # 例如：第一個時間點的結果
    """
    if all_results:
        first_time_key = list(all_results.keys())[0]  # the first time step
        first_result = all_results[first_time_key]
        df_integrated = first_result['df_integrated']
        interpolated_grid = first_result['interpolated_grid']  
        ds = first_result['dataset']
        
        print(f"\nExample: Accessing results for first time point ({first_time_key}):")
        print(f"    DataFrame shape: {df_integrated.shape}")
        print(f"    Grid shape: {interpolated_grid.shape}")
        print(f"    Dataset variables: {list(ds.data_vars.keys())}")
    """

    # 印出所有可用的變數供除錯使用
    """
    print(f"\nAvailable variables:")
    print(f"    all_results: {type(all_results)} with {len(all_results)} time points")
    print(f"    grid_info: {type(grid_info)} with keys {list(grid_info.keys())}")
    """
    
    # breakpoint()  # 最終結果檢查

#===========================================================================================