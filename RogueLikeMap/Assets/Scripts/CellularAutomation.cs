using UnityEngine;
using System.Collections;

public class CellularAutomation : MonoBehaviour {

    private const int width = 100;
    private const int height = 100;
    private int[,] map = new int[width, height];
    [Range(0, 1000)]
    public int roomsArea = 680;
    private void Update() {
        if (Input.GetMouseButtonDown(0)) {
            GenerateMap();
        }
    }

    private void RandomInit() {
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                if (i == 0 || i == width - 1 || j == 0 || j == height - 1) {
                    map[i, j] = 1;
                } else {
                    map[i, j] = Random.Range(0, 1000) < roomsArea ? 0 : 1;
                }
            }
        }
    }

    private int GetSurroundingWalls(int x, int y) {
        int count = 0;
        for (int i = x - 1; i <= x + 1; i++) {
            for (int j = y - 1; j <= y + 1; j++) {
                if (i >= 0 && i < width && j >= 0 && j < height) {
                    if (i != x || j != y) {
                        count += map[i, j];
                    }
                } else {
                    ++count;
                }
            }
        }
        return count;
    }

    private void GenerateMap() {
        RandomInit();
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                int wallCount = GetSurroundingWalls(i, j);
                if (wallCount > 4) {
                    map[i, j] = 1;
                }
            }
        }
    }

    private void OnDrawGizmos() {
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                Gizmos.color = map[i, j] == 1 ? Color.black : Color.white;
                Gizmos.DrawCube(new Vector2(i + 0.5f, j + 0.5f), Vector2.one);
            }
        }
    }
}