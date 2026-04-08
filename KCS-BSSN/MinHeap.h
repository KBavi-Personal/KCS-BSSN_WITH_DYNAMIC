// credit https://www.digitalocean.com/community/tutorials/min-heap-binary-tree
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "Triple.h"
class MinHeap {
public:
	std::vector<Triple<int, float, int>> arr; //value1, key, extra-value
	// Current Size of the Heap
	int size;
	// Maximum capacity of the heap
	int capacity;

	int parent(int i) { return (i - 1) / 2; }
	int left_child(int i) { return (2 * i + 1); }
	int right_child(int i) { return (2 * i + 2); }
	Triple<int, float, int> get_min() { return arr[0]; }

	void insert_minheap(Triple<int, float, int> element);
	void heapify(int index);
	void update_element(int index);
	void delete_minimum();
	void delete_element(int index);
	void print_heap();
	int find_element_index_value(int value);
	int find_element_index_value_extra_value(int value, int extra_value);
	MinHeap(int capacity);
	~MinHeap();

};
