package com.Innospectra.NanoScan;

import android.annotation.SuppressLint;
import android.content.Context;

import java.util.Arrays;

final class NoodleMoisturePredictor {
    static final class Result {
        final double moisture;
        final String status;

        Result(double moisture, String status) {
            this.moisture = moisture;
            this.status = status;
        }
    }

    private static final int SAVGOL_WINDOW = 11;
    private static final int SAVGOL_HALF_WINDOW = 5;
    private static final double EDGE_TOLERANCE_NM = 2.0;

    private final NoodleMoistureModel model;

    private NoodleMoisturePredictor(NoodleMoistureModel model) {
        this.model = model;
    }

    static NoodleMoisturePredictor load(Context context) throws Exception {
        return new NoodleMoisturePredictor(NoodleMoistureModel.load(context));
    }

    Result predict(double[] sourceWavelengths, double[] sourceAbsorbance) throws Exception {
        if (sourceWavelengths == null || sourceAbsorbance == null || sourceWavelengths.length != sourceAbsorbance.length) {
            throw new IllegalArgumentException("光譜資料長度不一致");
        }

        Spectrum spectrum = cleanAndSort(sourceWavelengths, sourceAbsorbance);
        if (spectrum.wavelengths.length < SAVGOL_WINDOW) {
            throw new IllegalArgumentException("有效波長點過少，請重新掃描");
        }

        double sourceMin = spectrum.wavelengths[0];
        double sourceMax = spectrum.wavelengths[spectrum.wavelengths.length - 1];
        double targetMin = model.wavelengths[0];
        double targetMax = model.wavelengths[model.wavelengths.length - 1];
        if (sourceMin - targetMin > EDGE_TOLERANCE_NM || targetMax - sourceMax > EDGE_TOLERANCE_NM) {
            throw new IllegalArgumentException("波長範圍不足，請重新掃描");
        }

        double[] resampled = pchipInterpolate(spectrum.wavelengths, spectrum.absorbance, model.wavelengths);
        double[] derivative = savitzkyGolayDerivative(resampled, 4.0);
        double[] normalized = snv(derivative);
        String status = signalStatus(spectrum.absorbance);
///計算公式
        double prediction = model.intercept;
        for (int i = 0; i < normalized.length; i++) {
            ///double scaled = (normalized[i] - model.xMean[i]) / model.xStd[i];
            double scaled = normalized[i] - model.xMean[i];
            prediction += scaled * model.coef[i];
        }
        return new Result(prediction, status);
    }

    @SuppressLint("DefaultLocale")
    private static String signalStatus(double[] absorbance) {
        double max = Double.NEGATIVE_INFINITY;
        for (double v : absorbance) {
            if (isFinite(v) && v > max) {
                max = v;
            }
        }
        if (max < 0.2) {
            return String.format("訊號過弱 Max=%.2f", max);
        }
        if (max > 3.0) {
            return String.format("訊號飽和 Max=%.2f", max);
        }
        return "正常";
    }

    private static double[] snv(double[] values) {
        double mean = 0.0;
        for (double value : values) {
            mean += value;
        }
        mean /= values.length;

        double variance = 0.0;
        for (double value : values) {
            double diff = value - mean;
            variance += diff * diff;
        }
        double std = Math.sqrt(variance / values.length);
        if (std == 0.0 || !Double.isFinite(std)) {
            throw new IllegalArgumentException("SNV 標準差為 0，無法正規化光譜");
        }

        double[] normalized = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            normalized[i] = (values[i] - mean) / std;
        }
        return normalized;
    }

    private static double[] savitzkyGolayDerivative(double[] y, double delta) {
        double[] result = new double[y.length];
        for (int i = 0; i < y.length; i++) {
            int start = i - SAVGOL_HALF_WINDOW;
            if (start < 0) {
                start = 0;
            }
            if (start > y.length - SAVGOL_WINDOW) {
                start = y.length - SAVGOL_WINDOW;
            }

            double s0 = SAVGOL_WINDOW;
            double s1 = 0.0;
            double s2 = 0.0;
            double s3 = 0.0;
            double s4 = 0.0;
            double t0 = 0.0;
            double t1 = 0.0;
            double t2 = 0.0;

            for (int j = 0; j < SAVGOL_WINDOW; j++) {
                double x = (start + j - i) * delta;
                double x2 = x * x;
                double value = y[start + j];
                s1 += x;
                s2 += x2;
                s3 += x2 * x;
                s4 += x2 * x2;
                t0 += value;
                t1 += value * x;
                t2 += value * x2;
            }

            double[][] matrix = {
                    {s0, s1, s2},
                    {s1, s2, s3},
                    {s2, s3, s4}
            };
            double[] rhs = {t0, t1, t2};
            result[i] = solve3x3(matrix, rhs)[1];
        }
        return result;
    }

    private static double[] solve3x3(double[][] a, double[] b) {
        double[][] m = new double[3][4];
        for (int r = 0; r < 3; r++) {
            System.arraycopy(a[r], 0, m[r], 0, 3);
            m[r][3] = b[r];
        }

        for (int p = 0; p < 3; p++) {
            int max = p;
            for (int r = p + 1; r < 3; r++) {
                if (Math.abs(m[r][p]) > Math.abs(m[max][p])) {
                    max = r;
                }
            }
            double[] tmp = m[p];
            m[p] = m[max];
            m[max] = tmp;

            double pivot = m[p][p];
            if (Math.abs(pivot) < 1e-12) {
                throw new IllegalArgumentException("Savitzky-Golay 計算失敗");
            }
            for (int c = p; c < 4; c++) {
                m[p][c] /= pivot;
            }
            for (int r = 0; r < 3; r++) {
                if (r == p) {
                    continue;
                }
                double factor = m[r][p];
                for (int c = p; c < 4; c++) {
                    m[r][c] -= factor * m[p][c];
                }
            }
        }
        return new double[]{m[0][3], m[1][3], m[2][3]};
    }

    private static double[] pchipInterpolate(double[] x, double[] y, double[] target) {
        int n = x.length;
        double[] h = new double[n - 1];
        double[] delta = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            h[i] = x[i + 1] - x[i];
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }

        double[] d = new double[n];
        if (n == 2) {
            d[0] = delta[0];
            d[1] = delta[0];
        } else {
            d[0] = pchipEndpointSlope(h[0], h[1], delta[0], delta[1]);
            d[n - 1] = pchipEndpointSlope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);
            for (int i = 1; i < n - 1; i++) {
                if (delta[i - 1] == 0.0 || delta[i] == 0.0 || Math.signum(delta[i - 1]) != Math.signum(delta[i])) {
                    d[i] = 0.0;
                } else {
                    double w1 = 2.0 * h[i] + h[i - 1];
                    double w2 = h[i] + 2.0 * h[i - 1];
                    d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
                }
            }
        }

        double[] out = new double[target.length];
        int interval = 0;
        for (int i = 0; i < target.length; i++) {
            double t = target[i];
            while (interval < n - 2 && t > x[interval + 1]) {
                interval++;
            }
            double hi = x[interval + 1] - x[interval];
            double s = (t - x[interval]) / hi;
            double s2 = s * s;
            double s3 = s2 * s;
            out[i] = (2.0 * s3 - 3.0 * s2 + 1.0) * y[interval]
                    + (s3 - 2.0 * s2 + s) * hi * d[interval]
                    + (-2.0 * s3 + 3.0 * s2) * y[interval + 1]
                    + (s3 - s2) * hi * d[interval + 1];
        }
        return out;
    }

    private static double pchipEndpointSlope(double h0, double h1, double m0, double m1) {
        double d = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
        if (Math.signum(d) != Math.signum(m0)) {
            return 0.0;
        }
        if (Math.signum(m0) != Math.signum(m1) && Math.abs(d) > Math.abs(3.0 * m0)) {
            return 3.0 * m0;
        }
        return d;
    }

    private static Spectrum cleanAndSort(double[] wavelengths, double[] absorbance) {
        Pair[] pairs = new Pair[wavelengths.length];
        int count = 0;
        for (int i = 0; i < wavelengths.length; i++) {
            if (isFinite(wavelengths[i]) && isFinite(absorbance[i])) {
                pairs[count++] = new Pair(wavelengths[i], absorbance[i]);
            }
        }
        pairs = Arrays.copyOf(pairs, count);
        Arrays.sort(pairs);

        double[] x = new double[count];
        double[] y = new double[count];
        int unique = 0;
        double previous = Double.NaN;
        for (int i = 0; i < count; i++) {
            if (i == 0 || pairs[i].x != previous) {
                x[unique] = pairs[i].x;
                y[unique] = pairs[i].y;
                previous = pairs[i].x;
                unique++;
            }
        }
        return new Spectrum(Arrays.copyOf(x, unique), Arrays.copyOf(y, unique));
    }

    private static final class Pair implements Comparable<Pair> {
        final double x;
        final double y;

        Pair(double x, double y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public int compareTo(Pair other) {
            return Double.compare(this.x, other.x);
        }
    }

    private static final class Spectrum {
        final double[] wavelengths;
        final double[] absorbance;

        Spectrum(double[] wavelengths, double[] absorbance) {
            this.wavelengths = wavelengths;
            this.absorbance = absorbance;
        }
    }

    private static boolean isFinite(double value) {
        return !Double.isNaN(value) && !Double.isInfinite(value);
    }
}
