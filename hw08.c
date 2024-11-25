#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>

#define MAX_V 1000  // Maximum number of vertices
#define MAX_N 100   // Maximum period

typedef struct {
    int vertex;
    int time_step;
    int weight;
} Node;

typedef struct {
    Node nodes[MAX_V * MAX_N];
    int size;
} PriorityQueue;

// Priority Queue functions
void initQueue(PriorityQueue* pq) {
    pq->size = 0;
}

void push(PriorityQueue* pq, Node node) {
    pq->nodes[pq->size++] = node;
    int i = pq->size - 1;
    while (i > 0 && pq->nodes[i].weight < pq->nodes[(i - 1) / 2].weight) {
        Node temp = pq->nodes[i];
        pq->nodes[i] = pq->nodes[(i - 1) / 2];
        pq->nodes[(i - 1) / 2] = temp;
        i = (i - 1) / 2;
    }
}

Node pop(PriorityQueue* pq) {
    Node minNode = pq->nodes[0];
    pq->nodes[0] = pq->nodes[--pq->size];
    int i = 0;
    while (2 * i + 1 < pq->size) {
        int minIndex = i;
        if (pq->nodes[2 * i + 1].weight < pq->nodes[minIndex].weight) {
            minIndex = 2 * i + 1;
        }
        if (2 * i + 2 < pq->size && pq->nodes[2 * i + 2].weight < pq->nodes[minIndex].weight) {
            minIndex = 2 * i + 2;
        }
        if (minIndex == i) break;
        Node temp = pq->nodes[i];
        pq->nodes[i] = pq->nodes[minIndex];
        pq->nodes[minIndex] = temp;
        i = minIndex;
    }
    return minNode;
}

bool isEmpty(PriorityQueue* pq) {
    return pq->size == 0;
}

// Graph and Dijkstra implementation
void dijkstra(int V, int N, int graph[MAX_V][MAX_V][MAX_N], int src, int dest) {
    int dist[MAX_V][MAX_N];
    for (int i = 0; i < V; i++)
        for (int j = 0; j < N; j++)
            dist[i][j] = INT_MAX;

    dist[src][0] = 0;
    PriorityQueue pq;
    initQueue(&pq);
    push(&pq, (Node){src, 0, 0});

    while (!isEmpty(&pq)) {
        Node current = pop(&pq);
        int u = current.vertex;
        int t = current.time_step;

        if (current.weight > dist[u][t]) continue;

        for (int v = 0; v < V; v++) {
            if (graph[u][v][t % N] != -1) {
                int next_t = (t + 1) % N;
                int new_weight = dist[u][t] + graph[u][v][t % N];
                if (new_weight < dist[v][next_t]) {
                    dist[v][next_t] = new_weight;
                    push(&pq, (Node){v, next_t, new_weight});
                }
            }
        }
    }

    // Find the shortest path to the destination at any time step
    int min_dist = INT_MAX;
    for (int t = 0; t < N; t++) {
        if (dist[dest][t] < min_dist) {
            min_dist = dist[dest][t];
        }
    }

    printf("Shortest distance from %d to %d: %d\n", src, dest, min_dist);
}

// Main function to parse input and run the algorithm
int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Usage: %s <graph_file>\n", argv[0]);
        return 1;
    }

    FILE* file = fopen(argv[1], "r");
    if (!file) {
        perror("Error opening file");
        return 1;
    }

    int V, N;
    fscanf(file, "%d %d", &V, &N);

    int graph[MAX_V][MAX_V][MAX_N];
    for (int i = 0; i < V; i++)
        for (int j = 0; j < V; j++)
            for (int k = 0; k < N; k++)
                graph[i][j][k] = -1;

    int vs, vt;
    while (fscanf(file, "%d %d", &vs, &vt) != EOF) {
        for (int i = 0; i < N; i++) {
            fscanf(file, "%d", &graph[vs][vt][i]);
        }
    }

    fclose(file);

    int src, dest;
    while (scanf("%d %d", &src, &dest) != EOF) {
        dijkstra(V, N, graph, src, dest);
    }

    return 0;
}
