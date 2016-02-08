package com.gazman.factor.matrix;

import com.gazman.factor.Logger;
import com.gazman.factor.VectorData;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Ilya Gazman on 2/8/2016.
 */
public class BitMatrix extends Logger {

    private BigInteger rows[];
    private BigInteger solutionRows[];
    private ArrayList<VectorData> vectorDatas;

    public ArrayList<ArrayList<VectorData>> solve(ArrayList<VectorData> vectorDatas) {
        this.vectorDatas = vectorDatas;
        rows = new BigInteger[vectorDatas.size()];
        solutionRows = new BigInteger[vectorDatas.size()];
        for (int i = 0; i < rows.length; i++) {
            rows[i] = new BigInteger(vectorDatas.get(i).toString(), 2);
            solutionRows[i] = BigInteger.ZERO.setBit(i);
        }

        log("init");
        printMatrix(rows);
        log("Init solution");
        printMatrix(solutionRows);

        HashMap<Integer, Object> map = new HashMap<>();
        for (int column = 0; column < rows.length; column++) {
            int selectedRow = -1;
            for (int row = 0; row < rows.length; row++) {
                if (!rows[row].testBit(column)) {
                    continue;
                }
                if (selectedRow == -1 && !map.containsKey(row)) {
                    selectedRow = row;
                    map.put(row, row);
                    continue;
                }
                if (selectedRow != -1) {
                    log("xor", selectedRow, row);
                    xor(row, selectedRow);
                    printMatrix(rows);
                    log("Solution");
                    printMatrix(solutionRows);
                }
            }
            for (int row = 0; row < selectedRow; row++) {
                if (!rows[row].testBit(column)) {
                    continue;
                }
                log("xor from bottom", selectedRow, row);
                xor(row, selectedRow);
                printMatrix(rows);
                log("Solution");
                printMatrix(solutionRows);
            }
        }

        log("Complete");
        printMatrix(rows);

        ArrayList<ArrayList<VectorData>> solutions = new ArrayList<>();
        for (int i = 0; i < rows.length; i++) {
            if (rows[i].equals(BigInteger.ZERO)) {
                solutions.add(createSolution(i));
            }
        }

        return solutions;
    }

    private ArrayList<VectorData> createSolution(int row) {
        ArrayList<VectorData> solution = new ArrayList<>();
        for (int column = 0; column < rows.length; column++) {
            if (solutionRows[row].testBit(column)) {
                solution.add(vectorDatas.get(column));
            }
        }
        return solution;
    }

    private void printMatrix(BigInteger rows[]) {
        if (!logsEnabled) {
            return;
        }
        for (int i = 0; i < rows.length; i++) {
            printRow(rows, i);
        }
        log();
    }

    private void printRow(BigInteger[] rows, int rowIndex) {
        if (!logsEnabled) {
            return;
        }
        BigInteger row = rows[rowIndex];
        for (int i = 0; i < rows.length; i++) {
            logInLine(row.testBit(i) ? 1 : 0);
        }
        log();
    }

    private void xor(int rowA, int rowB) {
        rows[rowA] = rows[rowA].xor(rows[rowB]);
        solutionRows[rowA] = solutionRows[rowA].xor(solutionRows[rowB]);
    }
}
