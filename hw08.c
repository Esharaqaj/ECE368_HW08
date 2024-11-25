#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <stdbool.h>

#define MAX_VERTICES 1000
#define MAX_EDGES 10000
#define INF INT_MAX

// Edge structure
typedef struct {
    int source, target, weights[10]; // Max period N is assumed to be <= 10
} Edge;

// Priority Queue Node
typedef struct {
    int vertex, time, dist;
} Node;

// Priority Queue
typedef struct {
    Node nodes[MAX_VERTICES * 10];
    int size;
} PriorityQueue;

// Graph structure
typedef struct {
    int V, N, edgeCount;
    Edge edges[MAX_EDGES];
} Graph;

// Function prototypes
void initPriorityQueue(PriorityQueue *pq);
void push(PriorityQueue *pq, int vertex, int time, int dist);
Node pop(PriorityQueue *pq);
bool isEmpty(PriorityQueue *pq);
void dijkstra(Graph *graph, int start, int end);
void readGraph(Graph *graph, const char *filename);

// Initialize Priority Queue
void initPriorityQueue(PriorityQueue *pq) {
    pq->size = 0;
}

// Push into Priority Queue
void push(PriorityQueue *pq, int vertex, int time, int dist) {
    pq->nodes[pq->size++] = (Node){vertex, time, dist};
    int i = pq->size - 1;

    // Bubble up
    while (i > 0 && pq->nodes[i].dist < pq->nodes[(i - 1) / 2].dist) {
        Node temp = pq->nodes[i];
        pq->nodes[i] = pq->nodes[(i - 1) / 2];
        pq->nodes[(i - 1) / 2] = temp;
        i = (i - 1) / 2;
    }
}

// Pop from Priority Queue
Node pop(PriorityQueue *pq) {
    Node root = pq->nodes[0];
    pq->nodes[0] = pq->nodes[--pq->size];
    int i = 0;

    // Bubble down
    while (2 * i + 1 < pq->size) {
        int smallest = i;
        int left = 2 * i + 1, right = 2 * i + 2;

        if (pq->nodes[left].dist < pq->nodes[smallest].dist) smallest = left;
        if (right < pq->size && pq->nodes[right].dist < pq->nodes[smallest].dist) smallest = right;

        if (smallest == i) break;

        Node temp = pq->nodes[i];
        pq->nodes[i] = pq->nodes[smallest];
        pq->nodes[smallest] = temp;
        i = smallest;
    }

    return root;
}

// Check if Priority Queue is empty
bool isEmpty(PriorityQueue *pq) {
    return pq->size == 0;
}

// Read Graph from file
void readGraph(Graph *graph, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%d %d", &graph->V, &graph->N);
    graph->edgeCount = 0;

    while (!feof(file)) {
        int source, target;
        fscanf(file, "%d %d", &source, &target);
        graph->edges[graph->edgeCount].source = source;
        graph->edges[graph->edgeCount].target = target;
        for (int i = 0; i < graph->N; i++) {
            fscanf(file, "%d", &graph->edges[graph->edgeCount].weights[i]);
        }
        graph->edgeCount++;
    }

    fclose(file);
}

// Dijkstra's Algorithm with Expanded Graph
void dijkstra(Graph *graph, int start, int end) {
    int dist[MAX_VERTICES][10]; // Distance array
    int prev[MAX_VERTICES][10]; // Previous node for path reconstruction
    bool visited[MAX_VERTICES][10] = {false}; // Visited array
    PriorityQueue pq;
    initPriorityQueue(&pq);

    // Initialize distances to infinity
    for (int i = 0; i < graph->V; i++) {
        for (int t = 0; t < graph->N; t++) {
            dist[i][t] = INF;
            prev[i][t] = -1;
        }
    }

    // Start from (start, 0)
    dist[start][0] = 0;
    push(&pq, start, 0, 0);

    while (!isEmpty(&pq)) {
        Node current = pop(&pq);
        int u = current.vertex;
        int t = current.time;

        if (visited[u][t]) continue;
        visited[u][t] = true;

        // Process neighbors
        for (int i = 0; i < graph->edgeCount; i++) {
            Edge *edge = &graph->edges[i];
            if (edge->source != u) continue;

            int v = edge->target;
            int nextTime = (t + 1) % graph->N;
            int weight = edge->weights[nextTime];
            int newDist = dist[u][t] + weight;

            if (newDist < dist[v][nextTime]) {
                dist[v][nextTime] = newDist;
                prev[v][nextTime] = u;
                push(&pq, v, nextTime, newDist);
            }
        }
    }

    // Find the minimum distance to the end node
    int minDist = INF, bestTime = -1;
    for (int t = 0; t < graph->N; t++) {
        if (dist[end][t] < minDist) {
            minDist = dist[end][t];
            bestTime = t;
        }
    }

    // Print the shortest path
    if (minDist == INF) {
        printf("No path found\n");
        return;
    }

    // Reconstruct path
    int path[MAX_VERTICES], pathSize = 0;
    for (int at = end, time = bestTime; at != -1; time = (time - 1 + graph->N) % graph->N) {
        path[pathSize++] = at;
        at = prev[at][time];
    }

    for (int i = pathSize - 1; i >= 0; i--) {
        printf("%d ", path[i]);
    }
    printf("\n");
}

// Main Function
int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <graph.txt>\n", argv[0]);
        return EXIT_FAILURE;
    }

    Graph graph;
    readGraph(&graph, argv[1]);

    int start, end;
    while (scanf("%d %d", &start, &end) == 2) {
        dijkstra(&graph, start, end);
    }

    return EXIT_SUCCESS;
}
