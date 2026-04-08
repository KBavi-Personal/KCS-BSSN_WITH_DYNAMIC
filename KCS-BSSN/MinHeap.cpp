// credit https://www.digitalocean.com/community/tutorials/min-heap-binary-tree
#include "MinHeap.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>

MinHeap::MinHeap(int _capacity) {
	capacity = _capacity;
	for (int i = 0; i < capacity; i++)
	{
		arr.push_back(Triple<int, float, int>(0, 0.0, 0));
	}
	size = 0;
}

void MinHeap::insert_minheap(Triple<int, float, int> element) {
	// Inserts an element to the min heap
	// We first add it to the bottom (last level)
	// of the tree, and keep swapping with it's parent
	// if it is lesser than it. We keep doing that until
	// we reach the root node. So, we will have inserted the
	// element in it's proper position to preserve the min heap property
	if (size == capacity) {
		fprintf(stderr, "Cannot insert %d. Heap is already full!\n", element.first);
		return;
	}
	// We can add it. Increase the size and add it to the end
	size++;
	arr[size - 1] = element;

	// Keep swapping until we reach the root
	int curr = size - 1;
	// As long as you aren't in the root node, and while the 
	// parent of the last element is greater than it
	while (curr > 0 && arr[parent(curr)].second > arr[curr].second) {
		// Swap
		Triple<int, float, int> temp = arr[parent(curr)];
		arr[parent(curr)] = arr[curr];
		arr[curr] = temp;
		// Update the current index of element
		curr = parent(curr);
	}
}

void MinHeap::heapify(int index) {
	// Rearranges the heap as to maintain
	// the min-heap property
	if (size <= 1)
		return;

	int left = left_child(index);
	int right = right_child(index);

	// Variable to get the smallest element of the subtree
	// of an element an index
	int smallest = index;

	// If the left child is smaller than this element, it is
	// the smallest
	if (left < size && arr[left].second < arr[index].second)
		smallest = left;

	// Similarly for the right, but we are updating the smallest element
	// so that it will definitely give the least element of the subtree
	if (right < size && arr[right].second < arr[smallest].second)
		smallest = right;

	// Now if the current element is not the smallest,
	// swap with the current element. The min heap property
	// is now satisfied for this subtree. We now need to
	// recursively keep doing this until we reach the root node,
	// the point at which there will be no change!
	if (smallest != index)
	{
		Triple<int, float, int> temp = arr[index];
		arr[index] = arr[smallest];
		arr[smallest] = temp;
		heapify(smallest);
	}
}

void MinHeap::update_element(int index)
{
	Triple<int, float, int> element = arr[index];
	delete_element(index);
	insert_minheap(element);

	//std::pair<int, float> first_element = arr[0];
	//arr[0] = arr[index];
	//arr[index] = first_element;
	////heapify(0)

}

void MinHeap::delete_minimum() {
	// Deletes the minimum element, at the root
	if (size == 0)
		return;

	Triple<int, float, int> last_element = arr[size - 1];

	// Update root value with the last element
	arr[0] = last_element;

	// Now remove the last element, by decreasing the size
	size--;
	// We need to call heapify(), to maintain the min-heap
	// property
	heapify(0);
}

void MinHeap::delete_element(int index) {
	Triple<int, float, int> min_element = get_min();
	// Deletes an element, indexed by index
	// Ensure that it's lesser than the current root
	arr[index] = Triple<int, float, int>(min_element.first, min_element.second - 1, min_element.third);

	// Now keep swapping, until we update the tree
	int curr = index;
	while (curr > 0 && arr[parent(curr)].second > arr[curr].second) {
		Triple<int, float, int> temp = arr[parent(curr)];
		arr[parent(curr)] = arr[curr];
		arr[curr] = temp;
		curr = parent(curr);
	}
	// Now simply delete the minimum element
	delete_minimum();
}

void MinHeap::print_heap() {
	// Simply print the array. This is an
	// inorder traversal of the tree
	printf("Min Heap:\n");
	for (int i = 0; i < size; i++) {
		printf("%d\t%f \n", arr[i].first, arr[i].second);
	}
	printf("\n");
}

MinHeap::~MinHeap() {
	if (size <= 0)
		return;
	arr.clear();
	size = 0;
	capacity = 0;

}
int MinHeap::find_element_index_value(int value) {
	if (size <= 0)
		return -1;
	for (int i = 0; i < size; i++)
		if (arr[i].first == value)
			return i;
	return -1;//not found
}
int MinHeap::find_element_index_value_extra_value(int value, int extra_value) {
	if (size <= 0)
		return -1;
	for (int i = 0; i < size; i++)
		if (arr[i].first == value && arr[i].third == extra_value)
			return i;
	return -1;//not found
}
