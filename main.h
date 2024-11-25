#ifndef MAIN_H
#define MAIN_H

typedef struct Edge {
    int target;
    int *weights;
    struct Edge *next;
} Edge;

typedef struct Graph {
    int V;
    int N;
    Edge **edges;
    int *edge_count;
} Graph;

Graph* parse_graph(const char *filename);
void dijkstra(Graph *graph, int start, int end);

#endif
