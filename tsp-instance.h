#ifndef TSP_INSTANCE_H
#define TSP_INSTANCE_H

typedef struct _TSPInstance TSPInstance;

TSPInstance *tspi_create(const char fileName[]);

int tspi_size(const TSPInstance *tspi);

int tspi_dist(const TSPInstance *tspi, int i, int j);

void tspi_free(TSPInstance *tspi);

#endif
