package com.gazman.factor.matrix;

import com.gazman.factor.Logger;
import com.gazman.factor.VectorData;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

/**
 * Created by Ilya Gazman on 2/8/2016.
 */
public class BitMatrix extends Logger {

    private BitSet rows[];
    private BitSet solutionRows[];
    private ArrayList<VectorData> vectorDatas;

    public ArrayList<ArrayList<VectorData>> solve(ArrayList<VectorData> vectorDatas) {
        log("Preparing to solve...");

        this.vectorDatas = vectorDatas;
        HashMap<Integer, Object> map = new HashMap<>();
        rows = new BitSet[vectorDatas.size()];
        solutionRows = new BitSet[vectorDatas.size()];
        for (int i = 0; i < rows.length; i++) {
            rows[i] = vectorDatas.get(i).vector;
            solutionRows[i] = new BitSet();
            solutionRows[i].set(i);
        }

        log("Solving...");

        for (int column = 0; column < rows.length; column++) {
            int selectedRow = -1;
            for (int row = 0; row < rows.length; row++) {
                if (!rows[row].get(column)) {
                    continue;
                }
                if (selectedRow == -1 && !map.containsKey(row)) {
                    selectedRow = row;
                    map.put(row, row);
                    continue;
                }
                if (selectedRow != -1) {
                    xor(row, selectedRow);
                }
            }
            for (int row = 0; row < selectedRow; row++) {
                if (!rows[row].get(column)) {
                    continue;
                }
                xor(row, selectedRow);
            }
        }

        log("Extracting solutions");

        ArrayList<ArrayList<VectorData>> solutions = new ArrayList<>();
        for (int i = 0; i < rows.length; i++) {
            if (rows[i].isEmpty()) {
                solutions.add(createSolution(i));
            }
        }

        log("Done");

        return solutions;
    }

    private ArrayList<VectorData> createSolution(int row) {
        ArrayList<VectorData> solution = new ArrayList<>();
        for (int column = 0; column < rows.length; column++) {
            if (solutionRows[row].get(column)) {
                solution.add(vectorDatas.get(column));
            }
        }
        return solution;
    }

    private void xor(int rowA, int rowB) {
        rows[rowA].xor(rows[rowB]);
        solutionRows[rowA].xor(solutionRows[rowB]);
    }
}
