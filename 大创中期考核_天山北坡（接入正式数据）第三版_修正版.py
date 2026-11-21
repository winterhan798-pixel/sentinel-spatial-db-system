"""真实 InSAR-地下水联合分析修正版。

相对于原始 Notebook：
- 只使用两类数据的共同观测窗口，不做边界外推；
- 对输入序列去线性趋势并做 Z-score 标准化；
- XWT 绘图提供 COI 遮罩、对数互功率和相位箭头；
- 批处理在 1 年周期上剔除 COI 边界，并使用单位复数向量圆平均；
- 自动跳过观测不足或方差近零的无效井。

原始 Notebook 文件保持不变。
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt
from scipy import signal
from scipy.interpolate import interp1d


WAVELET = "cmor1.5-1.0"
DT = 1.0 / 12.0
PERIODS_TARGET = np.logspace(np.log2(0.25), np.log2(2.5), 64, base=2)


def _zscore_detrend(values):
    detrended = signal.detrend(np.asarray(values, dtype=float))
    std = np.std(detrended)
    if not np.isfinite(std) or std < 1e-12:
        return None
    return (detrended - np.mean(detrended)) / std


def _common_monthly_series(sub_insar, sub_gw):
    """在共同观测窗口内插值，禁止使用 extrapolate。"""
    sub_insar = sub_insar.sort_values("Date").drop_duplicates("Date")
    sub_gw = sub_gw.sort_values("Date").drop_duplicates("Date")
    common_start = max(sub_insar["Date"].min(), sub_gw["Date"].min())
    common_end = min(sub_insar["Date"].max(), sub_gw["Date"].max())
    if common_end <= common_start:
        return None

    monthly_dates = pd.date_range(common_start, common_end, freq="ME")
    if len(monthly_dates) < 24:
        return None

    insar_days = (sub_insar["Date"] - common_start).dt.days.to_numpy()
    gw_days = (sub_gw["Date"] - common_start).dt.days.to_numpy()
    target_days = (monthly_dates - common_start).days.to_numpy()
    if target_days[0] < insar_days[0] or target_days[-1] > insar_days[-1]:
        return None
    if target_days[0] < gw_days[0] or target_days[-1] > gw_days[-1]:
        return None

    f_insar = interp1d(insar_days, sub_insar["InSAR_mm"].to_numpy(), kind="linear")
    f_gw = interp1d(gw_days, sub_gw["GW_level_m"].to_numpy(), kind="linear")
    return monthly_dates, f_insar(target_days), f_gw(target_days)


def _wavelet_xwt(insar_norm, gw_norm, dt=DT):
    center_frequency = pywt.central_frequency(WAVELET)
    scales = center_frequency / (PERIODS_TARGET * dt)
    cwt_insar, freqs = pywt.cwt(
        insar_norm, scales, WAVELET, sampling_period=dt
    )
    cwt_gw, _ = pywt.cwt(gw_norm, scales, WAVELET, sampling_period=dt)
    periods = 1.0 / freqs
    xwt = cwt_insar * np.conj(cwt_gw)
    power = np.abs(xwt)
    return periods, xwt, power


def plot_well_xwt(well_id, df_insar, df_gw, output_path=None):
    """绘制单井 XWT；返回图对象和中间结果。"""
    sub_insar = df_insar[df_insar["Well_ID"] == well_id]
    sub_gw = df_gw[df_gw["Well_ID"] == well_id]
    series = _common_monthly_series(sub_insar, sub_gw)
    if series is None:
        raise ValueError(f"{well_id}: 共同观测窗口不足或无法安全插值")

    dates, y_insar, y_gw = series
    insar_norm = _zscore_detrend(y_insar)
    gw_norm = _zscore_detrend(y_gw)
    if insar_norm is None or gw_norm is None:
        raise ValueError(f"{well_id}: 序列方差过小，无法进行小波分析")

    periods, xwt, power = _wavelet_xwt(insar_norm, gw_norm)
    log2_power = np.log2(power + 1e-8)
    time_years = dates.year + (dates.month - 0.5) / 12.0
    n = len(dates)
    edge_distance = np.minimum(np.arange(n), np.arange(n)[::-1]) * DT
    coi_period = edge_distance / np.sqrt(2.0)

    fig, ax = plt.subplots(figsize=(12, 6))
    levels = np.linspace(
        np.percentile(log2_power, 5), np.percentile(log2_power, 99), 20
    )
    contour = ax.contourf(
        time_years,
        periods,
        log2_power,
        levels=levels,
        cmap="jet",
        extend="both",
    )
    fig.colorbar(contour, ax=ax, label=r"交叉小波对数能量 $\log_2(|W_{xy}|)$")

    # 这是探索性分位数阈值，不冒充 AR(1) 统计意义上的 95% 检验。
    threshold = np.percentile(log2_power, 85)
    ax.contour(
        time_years,
        periods,
        log2_power,
        levels=[threshold],
        colors="black",
        linewidths=1.5,
    )
    ax.plot(time_years, coi_period, "k--", linewidth=1.3, label="COI 边界")
    ax.fill_between(
        time_years,
        coi_period,
        periods.max(),
        color="white",
        alpha=0.45,
    )

    t_idx = np.arange(2, n - 2, 3)
    p_idx = np.arange(4, len(periods) - 4, 5)
    t_grid, p_grid = np.meshgrid(time_years[t_idx], periods[p_idx])
    phase = np.angle(xwt[np.ix_(p_idx, t_idx)])
    mask = log2_power[np.ix_(p_idx, t_idx)] > np.percentile(log2_power, 65)
    ax.quiver(
        t_grid[mask],
        p_grid[mask],
        np.cos(phase)[mask],
        np.sin(phase)[mask],
        units="inches",
        scale=4.5,
        width=0.018,
        color="black",
        pivot="mid",
    )

    ax.axhline(1.0, color="white", linestyle=":", linewidth=1.6)
    ax.set_yscale("log", base=2)
    ax.set_ylim(periods.max(), periods.min())
    ax.set_xlabel("年份")
    ax.set_ylabel("周期（年）")
    ax.set_title(f"{well_id} 真实数据 XWT（去趋势标准化 + COI）")
    ax.legend(loc="upper right")
    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=180, bbox_inches="tight")
    return fig, {"dates": dates, "periods": periods, "xwt": xwt, "power": power}


def process_real_well(well_id, df_insar, df_gw, dt=DT):
    """处理单井并提取 1 年周期相位滞后。"""
    sub_insar = df_insar[df_insar["Well_ID"] == well_id].sort_values("Date")
    sub_gw = df_gw[df_gw["Well_ID"] == well_id].sort_values("Date")
    series = _common_monthly_series(sub_insar, sub_gw)
    if series is None:
        return None

    dates, y_insar, y_gw = series
    insar_norm = _zscore_detrend(y_insar)
    gw_norm = _zscore_detrend(y_gw)
    if insar_norm is None or gw_norm is None:
        return None

    periods, xwt, _ = _wavelet_xwt(insar_norm, gw_norm, dt)
    one_year_idx = int(np.argmin(np.abs(periods - 1.0)))
    xwt_1yr = xwt[one_year_idx, :]

    # 1 年周期的 e-folding COI 约为 sqrt(2) 年，按月剔除首尾边界。
    coi_margin = int(np.ceil(np.sqrt(2.0) / dt))
    if len(xwt_1yr) > 2 * coi_margin:
        valid_xwt = xwt_1yr[coi_margin:-coi_margin]
    else:
        return None
    unit_phasors = valid_xwt / (np.abs(valid_xwt) + 1e-12)
    avg_phase = np.angle(np.mean(unit_phasors))
    time_lag_days = (-avg_phase / (2.0 * np.pi)) * 365.0
    time_lag_days %= 365.0

    shift_steps = int(round(time_lag_days / (365.0 * dt)))
    if 0 < shift_steps < len(y_insar) - 12:
        aligned_insar = y_insar[shift_steps:]
        aligned_gw = y_gw[:-shift_steps]
    else:
        aligned_insar, aligned_gw = y_insar, y_gw

    ske = np.std(signal.detrend(aligned_insar)) / (
        np.std(signal.detrend(aligned_gw)) * 1000.0
    )
    t_rel = (dates - dates.iloc[0]).dt.days.to_numpy() / 365.0
    design = np.vstack(
        [t_rel, np.cos(2 * np.pi * t_rel), np.sin(2 * np.pi * t_rel), np.ones(len(t_rel))]
    ).T
    fit = np.linalg.lstsq(design, y_insar, rcond=None)[0]

    first = sub_insar.iloc[0]
    return {
        "Well_ID": well_id,
        "Well_Name": first.get("Well_Name", ""),
        "Aquifer_Type": first.get("Aquifer_Type", ""),
        "Lon": first.get("Lon", np.nan),
        "Lat": first.get("Lat", np.nan),
        "v_fit": fit[0],
        "A1_fit": np.hypot(fit[1], fit[2]),
        "Time_Lag_days": time_lag_days,
        "Ske_cal": ske,
        "n_months": len(dates),
        "common_start": dates.iloc[0].date().isoformat(),
        "common_end": dates.iloc[-1].date().isoformat(),
    }


def run(input_dir=".", output_csv="real_wells_summary_fixed.csv"):
    input_dir = Path(input_dir)
    df_insar = pd.read_csv(input_dir / "real_insar_ts.csv")
    df_gw = pd.read_csv(input_dir / "real_gw_ts.csv")
    df_insar["Date"] = pd.to_datetime(df_insar["Date"])
    df_gw["Date"] = pd.to_datetime(df_gw["Date"])

    common_wells = sorted(
        set(df_insar["Well_ID"].dropna()).intersection(df_gw["Well_ID"].dropna())
    )
    rows = [
        result
        for well_id in common_wells
        if (result := process_real_well(well_id, df_insar, df_gw)) is not None
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(input_dir / output_csv, index=False)
    print(f"有效井数: {len(summary)} / {len(common_wells)}")
    print(f"结果已写入: {input_dir / output_csv}")
    return summary, df_insar, df_gw


if __name__ == "__main__":
    summary, insar, gw = run()
    if not summary.empty:
        well_id = summary.iloc[0]["Well_ID"]
        plot_well_xwt(well_id, insar, gw, "sample_well_xwt_fixed.png")
        plt.show()
