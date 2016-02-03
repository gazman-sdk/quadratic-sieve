package com.gazman.factor;

import java.util.ArrayList;

/**
 * Created by Ilya Gazman on 1/30/2016.
 */
public class Matrix extends BaseFactor {
    private boolean[][] matrix;
    private boolean[][] solutionMatrix;

    public ArrayList<Object> solve(ArrayList<VectorData> vectorDatas) {
        matrix = new boolean[vectorDatas.size()][vectorDatas.get(0).vector.length];
        for (int i = 0; i < matrix.length; i++) {
            matrix[i] = vectorDatas.get(i).vector.clone();
        }
        int max = Math.max(matrix.length, matrix[0].length);
        solutionMatrix = new boolean[max][max];
        for (int i = 0; i < solutionMatrix.length; i++) {
            solutionMatrix[i][i] = true;
        }
        solve();

        ArrayList<Object> solutions = new ArrayList<>();

        for (int i = 0; i < matrix.length; i++) {
            boolean fail = false;
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j]) {
                    fail = true;
                    break;
                }
            }
            if (!fail) {
                solutions.add(solutionMatrix[i]);
            }
        }

        return solutions;
    }

    public void solve() {
        log("Solving");
        printMatrix(matrix);
        log("*************************");
        for (int i = 0; i < matrix.length; i++) {
            if (i < matrix[i].length && !matrix[i][i]) {
                swap(i);
            }
            if (i < matrix[i].length && matrix[i][i]) {
                clear(i);
            }
        }

        for (int i = matrix.length - 1; i > 0; i--) {
            if (i < matrix[i].length && matrix[i][i]) {
                clearUp(i);
            }
        }
    }

    private void printMatrix(boolean[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            logInLine(i + "> [");
            for (int j = 0; j < matrix[i].length; j++) {
                logInLine(matrix[i][j] ? 1 : 0);
            }
            logInLine("]");
            log();
        }
        log();
    }

    private void swap(int i) {
        for (int j = i + 1; j < matrix.length; j++) {
            if (matrix[j][i]) {
                swap(matrix, i, j);
                swap(solutionMatrix, i, j);
                break;
            }
        }
    }

    private void swap(boolean[][] matrix, int i, int j) {
        if(logsEnabled){
            log("Swap before", i, j, getMatrixName(matrix));
            printMatrix(matrix);
        }
        boolean[] tmp = matrix[j];
        matrix[j] = matrix[i];
        matrix[i] = tmp;
        if(logsEnabled){
            log("Swap after", i, j, getMatrixName(matrix));
            printMatrix(matrix);
        }
    }

    private String getMatrixName(boolean[][] matrix) {
        return matrix == this.matrix ? "Original" : "Solution";
    }

    private void clearUp(int i) {
        for (int j = i - 1; j >= 0; j--) {
            if (matrix[j][i]) {
                xor(matrix, i, j);
                xor(solutionMatrix, i, j);
            }
        }
    }

    private void clear(int i) {
        for (int j = i + 1; j < matrix.length; j++) {
            if (matrix[j][i]) {
                xor(matrix, i, j);
                xor(solutionMatrix, i, j);
            }
        }
    }

    private void xor(boolean[][] matrix, int i, int j) {
        if(logsEnabled){
            log("xor befor", i, j, getMatrixName(matrix));
            printMatrix(matrix);
        }
        for (int l = 0; l < matrix[j].length; l++) {
            matrix[j][l] = matrix[j][l] != matrix[i][l];
        }
        if(logsEnabled){
            log("xor after", i, j, getMatrixName(matrix));
            printMatrix(matrix);
        }
    }
}
