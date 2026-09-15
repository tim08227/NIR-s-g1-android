package com.Innospectra.NanoScan;

import android.content.Context;

import org.json.JSONArray;
import org.json.JSONObject;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;

final class NoodleMoistureModel {
    final double[] wavelengths;
    final double[] xMean;
    ///final double[] xStd;
    final double[] coef;
    final double intercept;

    private NoodleMoistureModel(
            double[] wavelengths,
            double[] xMean,
            ///double[] xStd,
            double[] coef,
            double intercept) {
        this.wavelengths = wavelengths;
        this.xMean = xMean;
        ///this.xStd = xStd;
        this.coef = coef;
        this.intercept = intercept;
    }

    static NoodleMoistureModel load(Context context) throws Exception {
        try (InputStream input = context.getAssets().open("noodle_moisture_model.json")) {
            ByteArrayOutputStream output = new ByteArrayOutputStream();
            byte[] buffer = new byte[4096];
            int read;
            while ((read = input.read(buffer)) != -1) {
                output.write(buffer, 0, read);
            }
            JSONObject json = new JSONObject(output.toString("UTF-8"));
            return new NoodleMoistureModel(
                    readDoubleArray(json.getJSONArray("wavelengths")),
                    readDoubleArray(json.getJSONArray("x_mean")),
                    ///readDoubleArray(json.getJSONArray("x_std")),
                    readDoubleArray(json.getJSONArray("coef")),
                    json.getDouble("intercept"));
        }
    }

    private static double[] readDoubleArray(JSONArray array) throws Exception {
        double[] values = new double[array.length()];
        for (int i = 0; i < array.length(); i++) {
            values[i] = array.getDouble(i);
        }
        return values;
    }
}
