#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>

#define INF INT_MAX

// Struct to represent an edge
typedef struct Edge {
    int target;
    int *weights;
} Edge;

// Struct to represent a graph
typedef struct Graph {
    int V;       // Number of vertices
    int N;       // Period of weights
    Edge **adj;  // Adjacency list
    int *adj_size; // Size of adjacency list for each vertex
} Graph;

// Priority queue node
typedef struct PriorityQueueNode {
    int vertex;
    int step;
    int weight;
} PriorityQueueNode;

// Min-heap for priority queue
typedef struct PriorityQueue {
    PriorityQueueNode *nodes;
    int size;
    int capacity;
} PriorityQueue;

// Function prototypes
Graph* create_graph(int V, int N);
void add_edge(Graph *graph, int src, int dest, int *weights);
void free_graph(Graph *graph);
PriorityQueue* create_priority_queue(int capacity);
void free_priority_queue(PriorityQueue *pq);
void push(PriorityQueue *pq, int vertex, int step, int weight);
PriorityQueueNode pop(PriorityQueue *pq);
int is_empty(PriorityQueue *pq);
int dijkstra(Graph *graph, int start, int end, int *path);
void handle_memory_error(void *ptr);
void free_all_resources(Graph *graph, PriorityQueue *pq);

// Main function
int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <graph_file>\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "r");
    if (!file) {
        perror("Error opening file");
        return 1;
    }

    // Parse the graph
    int V, N;
    if (fscanf(file, "%d %d", &V, &N) != 2 || V <= 0 || N <= 0) {
        fprintf(stderr, "Invalid graph header\n");
        fclose(file);
        return 1;
    }

    Graph *graph = create_graph(V, N);
    if (!graph) {
        fclose(file);
        return 1;
    }

    int src, dest;
    while (fscanf(file, "%d %d", &src, &dest) == 2) {
        int *weights = (int *)malloc(N * sizeof(int));
        handle_memory_error(weights);
        for (int i = 0; i < N; i++) {
            if (fscanf(file, "%d", &weights[i]) != 1) {
                fprintf(stderr, "Invalid edge weights\n");
                free(weights);
                free_graph(graph);
                fclose(file);
                return 1;
            }
        }
        add_edge(graph, src, dest, weights);
    }
    fclose(file);

    // Handle queries
    char query[256];
    while (fgets(query, sizeof(query), stdin)) {
        int start, end;
        if (sscanf(query, "%d %d", &start, &end) != 2) {
            fprintf(stderr, "Invalid query format\n");
            continue;
        }

        int *path = (int *)malloc(V * sizeof(int));
        handle_memory_error(path);
        int path_length = dijkstra(graph, start, end, path);
        if (path_length == -1) {
            printf("No path found\n");
        } else {
            for (int i = 0; i < path_length; i++) {
                printf("%d ", path[i]);
            }
            printf("\n");
        }
        free(path);
    }

    free_graph(graph);
    return 0;
}

// Create a graph
Graph* create_graph(int V, int N) {
    Graph *graph = (Graph *)malloc(sizeof(Graph));
    handle_memory_error(graph);
    graph->V = V;
    graph->N = N;
    graph->adj = (Edge **)malloc(V * sizeof(Edge *));
    handle_memory_error(graph->adj);
    graph->adj_size = (int *)calloc(V, sizeof(int));
    handle_memory_error(graph->adj_size);
    for (int i = 0; i < V; i++) {
        graph->adj[i] = NULL;
    }
    return graph;
}

// Add an edge to the graph
void add_edge(Graph *graph, int src, int dest, int *weights) {
    graph->adj[src] = (Edge *)realloc(graph->adj[src], (graph->adj_size[src] + 1) * sizeof(Edge));
    handle_memory_error(graph->adj[src]);
    graph->adj[src][graph->adj_size[src]].target = dest;
    graph->adj[src][graph->adj_size[src]].weights = weights;
    graph->adj_size[src]++;
}

// Free the graph
void free_graph(Graph *graph) {
    if (!graph) return;
    for (int i = 0; i < graph->V; i++) {
        for (int j = 0; j < graph->adj_size[i]; j++) {
            free(graph->adj[i][j].weights);
        }
        free(graph->adj[i]);
    }
    free(graph->adj);
    free(graph->adj_size);
    free(graph);
}

// Create a priority queue
PriorityQueue* create_priority_queue(int capacity) {
    PriorityQueue *pq = (PriorityQueue *)malloc(sizeof(PriorityQueue));
    handle_memory_error(pq);
    pq->nodes = (PriorityQueueNode *)malloc(capacity * sizeof(PriorityQueueNode));
    handle_memory_error(pq->nodes);
    pq->size = 0;
    pq->capacity = capacity;
    return pq;
}

// Free the priority queue
void free_priority_queue(PriorityQueue *pq) {
    if (pq) {
        free(pq->nodes);
        free(pq);
    }
}

// Push an element into the priority queue
void push(PriorityQueue *pq, int vertex, int step, int weight) {
    int i = pq->size++;
    while (i > 0 && pq->nodes[(i - 1) / 2].weight > weight) {
        pq->nodes[i] = pq->nodes[(i - 1) / 2];
        i = (i - 1) / 2;
    }
    pq->nodes[i] = (PriorityQueueNode){vertex, step, weight};
}

// Pop the minimum element from the priority queue
PriorityQueueNode pop(PriorityQueue *pq) {
    PriorityQueueNode min = pq->nodes[0];
    PriorityQueueNode last = pq->nodes[--pq->size];
    int i = 0;
    while (2 * i + 1 < pq->size) {
        int child = 2 * i + 1;
        if (child + 1 < pq->size && pq->nodes[child + 1].weight < pq->nodes[child].weight) {
            child++;
        }
        if (last.weight <= pq->nodes[child].weight) break;
        pq->nodes[i] = pq->nodes[child];
        i = child;
    }
    pq->nodes[i] = last;
    return min;
}

// Check if the priority queue is empty
int is_empty(PriorityQueue *pq) {
    return pq->size == 0;
}

// Handle memory allocation errors
void handle_memory_error(void *ptr) {
    if (!ptr) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
}

// Dijkstra's algorithm for periodic weights
int dijkstra(Graph *graph, int start, int end, int *path) {
    int V = graph->V, N = graph->N;

    // Dynamically allocate dist and prev arrays
    int **dist = (int **)malloc(V * sizeof(int *));
    int **prev = (int **)malloc(V * sizeof(int *));
    handle_memory_error(dist);
    handle_memory_error(prev);
    for (int i = 0; i < V; i++) {
        dist[i] = (int *)malloc(N * sizeof(int));
        prev[i] = (int *)malloc(N * sizeof(int));
        handle_memory_error(dist[i]);
        handle_memory_error(prev[i]);
        for (int j = 0; j < N; j++) {
            dist[i][j] = INF;
            prev[i][j] = -1;
        }
    }
    dist[start][0] = 0;

    PriorityQueue *pq = create_priority_queue(V * N);
    push(pq, start, 0, 0);

    while (!is_empty(pq)) {
        PriorityQueueNode node = pop(pq);
        int u = node.vertex;
        int step = node.step;

        if (u == end) {
            // Reconstruct path
            int length = 0;
            int v = u, s = step;
            while (v != -1) {
                path[length++] = v;
                int temp = prev[v][s];
                s = (s - 1 + N) % N;
                v = temp;
            }
            for (int i = 0; i < length / 2; i++) {
                int tmp = path[i];
                path[i] = path[length - 1 - i];
                path[length - 1 - i] = tmp;
            }
            free_priority_queue(pq);
            for (int i = 0; i < V; i++) {
                free(dist[i]);
                free(prev[i]);
            }
            free(dist);
            free(prev);
            return length;
        }

        for (int i = 0; i < graph->adj_size[u]; i++) {
            Edge edge = graph->adj[u][i];
            int v = edge.target;
            int next_step = (step + 1) % N;
            int weight = edge.weights[step];
            if (dist[u][step] + weight < dist[v][next_step]) {
                dist[v][next_step] = dist[u][step] + weight;
                prev[v][next_step] = u;
                push(pq, v, next_step, dist[v][next_step]);
            }
        }
    }

    free_priority_queue(pq);
    for (int i = 0; i < V; i++) {
        free(dist[i]);
        free(prev[i]);
    }
    free(dist);
    free(prev);
    return -1;
}
