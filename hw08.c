#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

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

// Struct to represent a priority queue node
typedef struct PriorityQueueNode {
    int node;
    int time;
    int weight;
} PriorityQueueNode;

// Struct to represent a priority queue
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
void resize_priority_queue(PriorityQueue *pq);
void free_priority_queue(PriorityQueue *pq);
void push(PriorityQueue *pq, int node, int time, int weight);
PriorityQueueNode pop(PriorityQueue *pq);
int is_empty(PriorityQueue *pq);
int dijkstra(Graph *graph, int start, int end, int *path);
void handle_memory_error(void *ptr);

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

// Resize the priority queue
void resize_priority_queue(PriorityQueue *pq) {
    pq->capacity *= 2;
    pq->nodes = (PriorityQueueNode *)realloc(pq->nodes, pq->capacity * sizeof(PriorityQueueNode));
    handle_memory_error(pq->nodes);
}

// Push an element into the priority queue
void push(PriorityQueue *pq, int node, int time, int weight) {
    if (pq->size == pq->capacity) {
        resize_priority_queue(pq);
    }

    int i = pq->size++;
    while (i > 0) {
        int parent = (i - 1) / 2;
        if (weight < pq->nodes[parent].weight ||
            (weight == pq->nodes[parent].weight && time < pq->nodes[parent].time)) {
            pq->nodes[i] = pq->nodes[parent];
            i = parent;
        } else {
            break;
        }
    }
    pq->nodes[i] = (PriorityQueueNode){node, time, weight};
}

// Pop the minimum element from the priority queue
PriorityQueueNode pop(PriorityQueue *pq) {
    PriorityQueueNode min = pq->nodes[0];
    PriorityQueueNode last = pq->nodes[--pq->size];
    int i = 0;
    while (2 * i + 1 < pq->size) {
        int left = 2 * i + 1, right = 2 * i + 2, smallest = left;

        if (right < pq->size &&
            (pq->nodes[right].weight < pq->nodes[left].weight ||
             (pq->nodes[right].weight == pq->nodes[left].weight &&
              pq->nodes[right].time < pq->nodes[left].time))) {
            smallest = right;
        }

        if (last.weight < pq->nodes[smallest].weight ||
            (last.weight == pq->nodes[smallest].weight && last.time <= pq->nodes[smallest].time)) {
            break;
        }

        pq->nodes[i] = pq->nodes[smallest];
        i = smallest;
    }
    pq->nodes[i] = last;
    return min;
}

// Check if the priority queue is empty
int is_empty(PriorityQueue *pq) {
    return pq->size == 0;
}

// Free the priority queue
void free_priority_queue(PriorityQueue *pq) {
    if (pq) {
        free(pq->nodes);
        free(pq);
    }
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

    // Initialize distance and parent arrays
    int **dist = (int **)malloc(V * sizeof(int *));
    int **parent = (int **)malloc(V * sizeof(int *));
    for (int i = 0; i < V; i++) {
        dist[i] = (int *)malloc(N * sizeof(int));
        parent[i] = (int *)malloc(N * sizeof(int));
        for (int j = 0; j < N; j++) {
            dist[i][j] = INF;
            parent[i][j] = -1;
        }
    }

    PriorityQueue *pq = create_priority_queue(V * N);
    push(pq, start, 0, 0);
    dist[start][0] = 0;

    while (!is_empty(pq)) {
        PriorityQueueNode curr = pop(pq);
        int u = curr.node;
        int t = curr.time;

        if (u == end) {
            // Reconstruct the path
            int length = 0;
            int node = u, time = t;
            while (node != -1) {
                path[length++] = node;
                int prev_time = (time - 1 + N) % N;
                int prev_node = parent[node][time];
                if (prev_node != -1 && graph->adj_size[prev_node] > 0) {
                    // Validate that the reconstructed edge is valid
                    int is_valid = 0;
                    for (int i = 0; i < graph->adj_size[prev_node]; i++) {
                        if (graph->adj[prev_node][i].target == node) {
                            is_valid = 1;
                            break;
                        }
                    }
                    if (!is_valid) {
                        fprintf(stderr, "Error: Invalid edge from %d to %d\n", prev_node, node);
                        return -1;
                    }
                }
                node = prev_node;
                time = prev_time;
            }
            for (int i = 0; i < length / 2; i++) {
                int temp = path[i];
                path[i] = path[length - i - 1];
                path[length - i - 1] = temp;
            }

            free_priority_queue(pq);
            for (int i = 0; i < V; i++) {
                free(dist[i]);
                free(parent[i]);
            }
            free(dist);
            free(parent);

            return length;
        }

        for (int i = 0; i < graph->adj_size[u]; i++) {
            Edge edge = graph->adj[u][i];
            int v = edge.target;
            int next_time = (t + 1) % N;
            int weight = edge.weights[t];
            if (dist[u][t] + weight < dist[v][next_time]) {
                dist[v][next_time] = dist[u][t] + weight;
                parent[v][next_time] = u;  // Correctly track the parent
                push(pq, v, next_time, dist[v][next_time]);
            }
        }
    }

    // If we reach here, there is no path to the destination
    free_priority_queue(pq);
    for (int i = 0; i < V; i++) {
        free(dist[i]);
        free(parent[i]);
    }
    free(dist);
    free(parent);

    return -1;
}