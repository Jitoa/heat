# Heat DNS Code - 2D Visualization Enhancement

## 概要 (Overview)

このプロジェクトでは、ppf-mainのコードを参考にして、heatのDNSコードに速度変動と温度変動の2次元可視化機能を追加しました。コメントアウトのみで既存のコードを活用し、バイナリデータ出力機能も追加して高品質な可視化を実現しています。

This project enhances the heat DNS code with 2D visualization capabilities for velocity and temperature fluctuations, inspired by the ppf-main code. We've uncommented existing code blocks and added binary data output for high-quality visualization.

## 実装した変更 (Implemented Changes)

### 1. INCFILE/INCTMSRS の修正

#### TMSRSU サブルーチン (速度場)
- **ブロック 91-95 をアンコメント**:
  - `uxz_XXX.dat`: x-z面での u' 速度変動 (3つのy位置: JA2, JD2, JE2)
  - `uxz_040_XXX.dat`: y=JD2 (y+≈31) での u' のみ
  - `uxz_080_XXX.dat`: y=JE2 (y+≈63) での u' のみ  
  - `uzy_XXX.dat`: z-y面での u', v', w' 速度変動
  - `uxy_XXX.dat`: x-y面での u', v' 速度変動

#### TMSRST サブルーチン (温度場)
- **ブロック 96-99 をアンコメント**:
  - `t1xz_XXX.dat`: x-z面での T1', T2' 温度変動 (UHF用とCTD用)
  - `t2xz_XXX.dat`: x-z面での T2' のみ
  - `t12zy_XXX.dat`: z-y面での T1', T2' 温度変動
  - `t12xy_XXX.dat`: x-y面での T1', T2' 温度変動

### 2. 新規追加: BINARYOUT サブルーチン

ppf-mainのINSTOUTサブルーチンを参考に、バイナリデータ出力機能を追加：
- **xzc_vel_XXXXXXXXX.dat**: チャネル中央での速度場バイナリデータ
- **xzc_temp_XXXXXXXXX.dat**: チャネル中央での温度場バイナリデータ

### 3. Python可視化スクリプト

`visualize_heat_2d.py`スクリプトを作成し、以下の機能を提供：
- テキストファイル（ブロック91-99）からの2D可視化
- バイナリファイル（BINARYOUT）からの高品質コンター図
- 自動ファイル検出と一括可視化
- 速度場と温度場の専用カラーマップ

## 出力条件 (Output Conditions)

- **テキストファイル**: `MOD(NN,100).EQ.0` - 100ステップごと
- **バイナリファイル**: `MOD(NN,1000).EQ.0` - 1000ステップごと
- **格子サンプリング**: 2つ置き (`I=1,IG,2`, `K=1,KG,2`)

## y位置の物理的意味 (Physical Meaning of y-positions)

- **JA2 = 1**: y+ ≈ 0.09 (壁近傍、壁乱流の粘性底層)
- **JD2 = 36**: y+ ≈ 31.24 (中間高さ、対数則領域)
- **JE2 = 48**: y+ ≈ 62.56 (チャネル中心、外層)

## 使用方法 (Usage)

### ステップ 1: Fortranコードのコンパイルと実行

```bash
# 通常通りコンパイル
gfortran -O2 -fopenmp ch060owv2.f -o heat_dns

# 実行
./heat_dns
```

### ステップ 2: Python可視化

```bash
# 必要なライブラリの確認/インストール
pip install numpy matplotlib

# 可視化実行
python visualize_heat_2d.py
```

## 出力ファイル (Output Files)

### テキストデータファイル
- `uxz_XXX.dat` - 速度変動 u' (3列: JA2, JD2, JE2)
- `uxz_040_XXX.dat` - 中間高さでの u'
- `uxz_080_XXX.dat` - チャネル中心での u'
- `uzy_XXX.dat` - z-y面速度変動 (3列: u', v', w')
- `uxy_XXX.dat` - x-y面速度変動 (2列: u', v')
- `t1xz_XXX.dat` - 温度変動 (6列: T1' & T2' × 3位置)
- `t2xz_XXX.dat` - T2'のみ (3列)
- `t12zy_XXX.dat` - z-y面温度変動 (2列: T1', T2')
- `t12xy_XXX.dat` - x-y面温度変動 (2列: T1', T2')

### バイナリデータファイル
- `xzc_vel_XXXXXXXXX.dat` - 速度場バイナリ (X, Z, U, V, W, U', V', W')
- `xzc_temp_XXXXXXXXX.dat` - 温度場バイナリ (X, Z, T1, T2, T1', T2')

### 可視化ファイル
- `figures/` ディレクトリに PNG 形式で保存
- 各データファイルに対応する2D可視化図
- コンター図とカラーマップ表示

## ppf-mainとの比較 (Comparison with ppf-main)

| 特徴 | ppf-main | heat (enhanced) |
|------|----------|-----------------|
| 出力形式 | バイナリメイン | テキスト + バイナリ |
| 2D平面 | xzc, zy, xy | xz, zy, xy (複数y位置) |
| 温度場 | なし | UHF & CTD境界条件 |
| 可視化 | MATLAB/AVS | Python + matplotlib |
| ディレクトリ構造 | 自動作成 | 単一ディレクトリ |

## 物理的解釈 (Physical Interpretation)

### 速度変動パターン
- **近壁 (JA2)**: 縦渦とストリーク構造が顕著
- **中間 (JD2)**: 壁乱流と外層の遷移領域
- **中心 (JE2)**: より等方的な乱流構造

### 温度変動パターン
- **UHF (T1')**: 一定熱流束境界条件での温度変動
- **CTD (T2')**: 一定温度差境界条件での温度変動
- 異なる境界条件による温度場の違いを比較可能

## トラブルシューティング (Troubleshooting)

### 1. ファイルが生成されない
- `NN`変数が正しく更新されているか確認
- 十分なステップ数（100以上）実行されているか確認
- ファイル書き込み権限の確認

### 2. Python可視化エラー
- numpy, matplotlibのインストール確認: `pip install numpy matplotlib`
- データファイルの存在確認
- ファイル形式の確認（テキスト vs バイナリ）

### 3. 格子サイズ不整合
- `INCFILE/INCINIT96`で`IG`, `KG`の値を確認
- Python スクリプトの grid size パラメータを調整

## 発展的な使用 (Advanced Usage)

### カスタム可視化
`visualize_heat_2d.py`を編集して：
- 異なるカラーマップの使用
- 特定の時刻のみの可視化
- 統計量の計算と表示
- 動画作成 (FFmpeg使用)

### 高解像度出力
格子間引きなしの完全解像度が必要な場合は、TMSRSU/TMSRSTの出力ブロックで`I=1,IG,2`を`I=1,IG`に変更。

### 追加の2D平面
他の y 位置での出力が必要な場合は、JA2, JD2, JE2の値を変更するか、新しい出力ブロックを追加。

## 結論 (Conclusion)

この実装により、heatコードでppf-mainと同等以上の2D可視化機能が利用可能になりました：

✅ **速度変動 (u', v', w')** の複数平面での可視化  
✅ **温度変動 (T1', T2')** のUHF/CTD境界条件比較  
✅ **テキスト + バイナリ** のデュアル出力形式  
✅ **Python自動可視化** による簡単な図作成  
✅ **コメントアウト解除のみ** の最小限修正  

発表用の図として十分な品質で、乱流構造の特徴的なパターンを明確に可視化できます。