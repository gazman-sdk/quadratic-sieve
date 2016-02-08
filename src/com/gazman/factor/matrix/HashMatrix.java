package com.gazman.factor.matrix;

import com.gazman.factor.Logger;
import com.gazman.factor.VectorData;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Created by Ilya Gazman on 1/30/2016.
 */
public class HashMatrix extends Logger {
    private static final Object CELL = new Object();
    private ArrayList<ConcurrentHashMap<Integer, Object>> columns = new ArrayList<>();
    private ArrayList<ConcurrentHashMap<Integer, Object>> rows = new ArrayList<>();
    private HashMatrix solutionMatrix;
    private ArrayList<VectorData> vectorDatas;

    public ArrayList<ArrayList<VectorData>> solve(ArrayList<VectorData> vectorDatas) {
        this.vectorDatas = vectorDatas;
        init(vectorDatas);
        solutionMatrix = new HashMatrix();
        solutionMatrix.setLogsEnabled(false);
        HashMap<Integer, Object> map = new HashMap<>();
        solutionMatrix.initWith1(Math.max(vectorDatas.size(), vectorDatas.get(0).vector.length));

        forceLog("Starting To Solve");

        for (int i = 0; i < columns.size(); i++) {
            ConcurrentHashMap<Integer, Object> mainColumn = columns.get(i);
            switch (mainColumn.size()){
                case 0:
                    continue;
                case 1:
                    map.put(mainColumn.keySet().iterator().next(), CELL);
                    continue;
            }

            Integer mainRow = null;

            for (Integer row :  mainColumn.keySet()) {
                if (!map.containsKey(row)) {
                    map.put(row, CELL);
                    mainRow = row;
                    break;
                }
            }

            if(mainRow != null){
                for (Integer row :  mainColumn.keySet()) {
                    if(row.equals(mainRow)){
                        continue;
                    }
                    for (Integer column : rows.get(mainRow).keySet()) {
                        xor(row,  column);
                    }
                    solutionMatrix.xorRows(mainRow,  row);
                }
            }
        }

        forceLog("Building solutions");

        ArrayList<ArrayList<VectorData>> solutions = new ArrayList<>();
        for (int i = 0; i < rows.size(); i++) {
            ConcurrentHashMap<Integer, Object> row = rows.get(i);
            if (row.size() == 0) {
                solutions.add(createSolution(i));
            }
        }

        forceLog("Solutions ready");

        return solutions;
    }

    private void xorRows(Integer mainRow, Integer row) {
        for (Integer localColumn : rows.get(mainRow).keySet()) {
            if(localColumn != null){
                xor(row,  localColumn);
            }
        }
    }

    private ArrayList<VectorData> createSolution(int row) {
        ArrayList<VectorData> solution = new ArrayList<>();
        for (Integer cell : solutionMatrix.rows.get(row).keySet()) {
            solution.add(vectorDatas.get(cell));
        }
        return solution;
    }

    private void initWith1(int size) {
        for (int i = 0; i < size; i++) {
            ConcurrentHashMap<Integer, Object> row = new ConcurrentHashMap<>();
            ConcurrentHashMap<Integer, Object> column = new ConcurrentHashMap<>();
            row.put(i, CELL);
            column.put(i, CELL);
            rows.add(row);
            columns.add(column);
        }
    }

    private void init(ArrayList<VectorData> vectorDatas) {
        for (int i = 0; i < vectorDatas.size(); i++) {
            rows.add(new ConcurrentHashMap<Integer, Object>());
        }
        int columnsCount = vectorDatas.get(0).vector.length;
        for (int i = 0; i < columnsCount; i++) {
            columns.add(new ConcurrentHashMap<Integer, Object>());
        }

        for (int row = 0; row < vectorDatas.size(); row++) {
            boolean[] cells = vectorDatas.get(row).vector;
            for (int cell = 0; cell < cells.length; cell++) {
                if (cells[cell]) {
                    add(cell, row);
                }
            }
        }
    }

    private void add(int row, int column) {
        rows.get(column).put(row, CELL);
        columns.get(row).put(column, CELL);
    }

    private void xor(int row, int column) {
        xor(rows.get(row), column);
        xor(columns.get(column), row);
    }

    private void xor(ConcurrentHashMap<Integer, Object> map, Integer cell) {
        if(map.containsKey(cell)){
            map.remove(cell);
        }
        else{
            map.put(cell, CELL);
        }
    }
}
