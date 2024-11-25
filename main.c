#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "main.h"

// Node for the priority queue
typedef struct Node {
    int vertex;
    int time;
    int weight;
} Node;

// Priority queue structure
typedef struct PriorityQueue {
    Node *nodes;
    int size;
    int capacity;
} PriorityQueue;

// Function to initialize the priority queue
PriorityQueue* create_priority_queue(int capacity) {
    PriorityQueue *pq = (PriorityQueue *)malloc(sizeof(PriorityQueue));
    pq->nodes = (Node *)malloc(capacity * sizeof(Node));
    pq->size = 0;
    pq->capacity = capacity;
    return pq;
}

// Function to add a node to the priority queue
void enqueue(PriorityQueue *pq, int vertex, int time, int weight) {
    pq->nodes[pq->size++] = (Node){vertex, time, weight};
    int i = pq->size - 1;
    while (i > 0 && pq->nodes[i].weight < pq->nodes[(i - 1) / 2].weight) {
        Node temp = pq->nodes[i];
        pq->nodes[i] = pq->nodes[(i - 1) / 2];
        pq->nodes[(i - 1) / 2] = temp;
        i = (i - 1) / 2;
    }
}

// Function to remove the node with the smallest weight from the queue
Node dequeue(PriorityQueue *pq) {
    Node root = pq->nodes[0];
    pq->nodes[0] = pq->nodes[--pq->size];
    int i = 0;
    while (i * 2 + 1 < pq->size) {
        int smallest = i;
        if (pq->nodes[i * 2 + 1].weight < pq->nodes[smallest].weight) smallest = i * 2 + 1;
        if (i * 2 + 2 < pq->size && pq->nodes[i * 2 + 2].weight < pq->nodes[smallest].weight) smallest = i * 2 + 2;
        if (smallest == i) break;
        Node temp = pq->nodes[i];
        pq->nodes[i] = pq->nodes[smallest];
        pq->nodes[smallest] = temp;
        i = smallest;
    }
    return root;
}

// Function to parse the graph input file
Graph* parse_graph(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Unable to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    int V, N;
    fscanf(file, "%d %d", &V, &N);

    Graph *graph = (Graph *)malloc(sizeof(Graph));
    graph->V = V;
    graph->N = N;
    graph->edges = (Edge **)malloc(V * sizeof(Edge *));
    graph->edge_count = (int *)calloc(V, sizeof(int));

    for (int i = 0; i < V; i++) {
        graph->edges[i] = NULL;
    }

    int vs, vt;
    while (fscanf(file, "%d %d", &vs, &vt) != EOF) {
        Edge *edge = (Edge *)malloc(sizeof(Edge));
        edge->target = vt;
        edge->weights = (int *)malloc(N * sizeof(int));
        for (int i = 0; i < N; i++) {
            fscanf(file, "%d", &edge->weights[i]);
        }
        edge->next = graph->edges[vs];
        graph->edges[vs] = edge;
        graph->edge_count[vs]++;
    }

    fclose(file);
    return graph;
}

// Dijkstra's algorithm for the expanded graph
void dijkstra(Graph *graph, int start, int end) {
    int V = graph->V;
    int N = graph->N;
    int **dist = (int **)malloc(V * sizeof(int *));
    int **prev = (int **)malloc(V * sizeof(int *));
    for (int i = 0; i < V; i++) {
        dist[i] = (int *)malloc(N * sizeof(int));
        prev[i] = (int *)malloc(N * sizeof(int));
        for (int j = 0; j < N; j++) {
            dist[i][j] = INT_MAX;
            prev[i][j] = -1;
        }
    }

    PriorityQueue *pq = create_priority_queue(V * N);
    enqueue(pq, start, 0, 0);
    dist[start][0] = 0;

    while (pq->size > 0) {
        Node node = dequeue(pq);
        int u = node.vertex;
        int t = node.time;

        if (u == end) break;

        for (Edge *edge = graph->edges[u]; edge != NULL; edge = edge->next) {
            int v = edge->target;
            int weight = edge->weights[t % N];
            int new_time = (t + 1) % N;
            int new_dist = dist[u][t] + weight;

            if (new_dist < dist[v][new_time]) {
                dist[v][new_time] = new_dist;
                prev[v][new_time] = u;
                enqueue(pq, v, new_time, new_dist);
            }
        }
    }

    // Output the path
    int t = 0, current = end;
    for (int i = 0; i < N; i++) {
        if (dist[end][i] < dist[end][t]) t = i;
    }

    int path[V], path_len = 0;
    while (current != -1) {
        path[path_len++] = current;
        current = prev[current][t];
        t = (t - 1 + N) % N;
    }

    for (int i = path_len - 1; i >= 0; i--) {
        printf("%d%c", path[i], i > 0 ? ' ' : '\n');
    }

    // Free resources
    free(pq->nodes);
    free(pq);
    for (int i = 0; i < V; i++) {
        free(dist[i]);
        free(prev[i]);
    }
    free(dist);
    free(prev);
}

// Main function
int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <graph file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    Graph *graph = parse_graph(argv[1]);

    int start, end;
    while (scanf("%d %d", &start, &end) == 2) {
        dijkstra(graph, start, end);
    }

    return EXIT_SUCCESS;
}
