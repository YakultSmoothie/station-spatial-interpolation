#  空間內插工具 (測站資料 -> 指定網格)

- 氣象測站資料空間內插Python工具，支援多種內插方法與彈性網格配置。
- create_station_spatial_interpolation.py [v1.0]

## 專案概述

本工具解決傳統網格內插程式的限制，提供以下功能：
- 彈性的網格位置與解析度配置
- 支援多種氣象變數（不僅限於降水資料）
- 多種內插方法與效能最佳化
- 多時間點批次處理
- UTC時間轉換功能
- 預覽視覺化輸出

## 核心功能

### 主要特色
- **彈性網格系統**：使用者自訂經緯度範圍與解析度
- **多變數支援**：支援任何具有經緯度與數值結構的氣象變數
- **時間序列處理**：多時間點批次處理
- **UTC轉換**：自動本地時間轉UTC時間

### 內插方法
- **1.1 Linear**：使用 scipy.griddata 的標準線性內插
- **1.2 Cubic**：三次樣條內插(Cubic spline interpolation)，產生平滑內插結果
- **1.3 Nearest**：最近鄰內插，適用於類別資料
- **2 IDW (距離反比權重法)**：基於距離的權重方法，可設定參數截段距離與權重次方。本程序直接計算，非使用import。
- **3 Kriging**：使用 pykrige.ok.OrdinaryKriging 的球面變異函數模型的克利金法

### 輸入
- 1. 測站基本資訊檔案(.txt格式)
- 2. 觀測資料檔案(.txt格式)
  
### 輸出
- **CSV**：內插使用的測站資料
- **NetCDF**：內插網格資料（CF相容，支援GrADS）
- **PNG**：測站觀測與內插網格的比較視覺化

## 安裝需求

### 必要套件
```bash
# 核心科學計算套件
pip install --user pandas numpy xarray scipy

# 安裝工作包 netCDF C 函式庫
pip3 install --user netCDF4

# 地理空間與視覺化
pip install --user matplotlib cartopy

# 選用：克利金內插法
pip install --user pykrige
```

### 開發環境（環境參考）
- Python 3.10.12
- OS: Ubuntu 22.04.3 LTS

## 使用方法

### 基本用法
```bash
python3 create_station_spatial_interpolation.py \
    -i ./data_station.txt \
    -i2 ./data_20250730_pp01.txt \
    -V PP01 \
    -T 2025073001,2025073002,2025073003
```

## 輸入資料格式

### 測站基本資訊檔案
```
station_id longitude latitude elevation
AAAAAA     120.7811  25.0711  756.0
BBBBBB     121.4419  24.9975  10.0
...
```

### 觀測資料檔案
```
station_id yyyymmddhh TX01 PP01
AAAAAA 2025073001 28.9 0.0
AAAAAA 2025073002 29.0 0.0
...
BBBBBB 2025073001 27.9 10.0
BBBBBB 2025073002 28.0 20.0
...
```

## 輸出結構

```
output/
├── csv/
│   ├── integrated_station_data_PP01_2025073001.csv
│   └── integrated_station_data_PP01_... ....png
├── nc/
│   ├── interpolated_PP01_2025073001.nc
│   └── interpolated_PP01_... ....nc
└── png/
    ├── interpolation_visualization_PP01_2025073001_linear.png
    └── interpolation_visualization_PP01_... ..._linear.png
```

## 命令列參數

### 必要參數
- `-i, --input`：測站基本資訊檔案路徑
- `-i2, --input2`：觀測資料檔案路徑

### 選用參數
| 參數 | 預設值 | 說明 |
|------|--------|------|
| `-V, --variable` | PP01 | 目標變數名稱 (變數名需對應觀測資料檔案) |
| `-T, --time` | 2025073001,2025073002 | 目標時間（逗號分隔） |
| `-LL, --lonlat` | 120.0,122.12,21.9,25.4 | 網格邊界（經度最小,經度最大,緯度最小,緯度最大） |
| `-n, --output_dir` | ./output | 輸出目錄路徑 |
| `-r, --resolution` | 0.03 | 網格解析度（度） |
| `-m, --method` | linear | 內插方法 |
| `--idw_power` | 3.0 | IDW冪次參數 |
| `--max_distance` | 50.0 | IDW最大搜尋距離（公里） |
| `-nd, --no_vis` | False | 停用視覺化輸出 |
| `-utc, --utc_offset` | None | UTC轉換偏移量（+N或-N小時）(轉輸出檔名為_UTC) |

## GrADS整合

NetCDF檔案符合CF標準並與GrADS相容。時間序列使用template選項：

```grads
dset ^./output/interpolated_PP01_%y4%m2%d2%h2.nc
options template
dtype netcdf
undef NaN
xdef 72 linear 120.0 0.03
ydef 114 linear 21.9 0.03  
zdef 1 linear 1 1
tdef 24 linear 01Z30jul2025 1hr
vars 1
    PP01=>r1  1  y,x  [mm/h]
endvars
```

## 技術細節

### 品質控制
- 測站-觀測資料匹配驗證

### 時間處理
- 支援24小時格式，hh為24時轉換至次日00:00
- 不同時區的彈性UTC轉換

## 使用範例

### 比較不同內插方法
```bash
# 線性內插法
python3 create_station_spatial_interpolation.py -i data_station.txt -i2 data_20250730_pp01.txt -V PP01 -T 2025073001,2025073002,2025073003 -m linear -n output_linear

# IDW內插法
python3 create_station_spatial_interpolation.py -i data_station.txt -i2 data_20250730_pp01.txt -V PP01 -T 2025073001,2025073002,2025073003 -m idw --idw_power 3 --max_distance 30 -n output_idw
```

### 後處理範例
```bash
# 後處理範例 - 3小時累積降水
grads -bcl "sample_plot.gs"
```

## 注意事項

### 遺失值
- 注意原始檔案是否有遺失值，在進入本程式前需要先預處理遺失值為NaN。
- 不同來源的資料使用的不同遺失值。

### 測站基本資訊
- 使用時需確保測站資訊檔案為最新且完整， 因為中央氣象署自動測站網隨時間變化，有新增/退役測站。
- 確認基本資訊與觀測資料間的測站代碼匹配

### 時區處理
- 測站資料通常使用當地時間（例如中央氣象署自動站使用UTC+8），與其他的結果比較時要注意對應的時間。

### 網格解析度
- 大型區域建議先使用較粗解析度測試
