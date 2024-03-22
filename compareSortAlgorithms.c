#include <time.h>
#include <stdlib.h>
#include <stdio.h>

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
	extraMemoryAllocated += sz;
	size_t* ret = malloc(sizeof(size_t) + sz);
	*ret = sz;
	printf("Extra memory allocated, size: %ld\n", sz);
	return &ret[1];
}

void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	printf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}

// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated
void heapify(int arr[], int n, int i) {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    if (l < n && arr[l] > arr[largest])
        largest = l;

    if (r < n && arr[r] > arr[largest])
        largest = r;

    if (largest != i) {
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;

        heapify(arr, n, largest);
    }
}

// main function to do heap sort
void heapSort(int arr[], int n) {
    int *heap = (int *)Alloc(sizeof(int) * n);

    memcpy(heap, arr, sizeof(int) * n);

    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(heap, n, i);

    for (int i = n - 1; i > 0; i--) {

        int temp = heap[0];
        heap[0] = heap[i];
        heap[i] = temp;

        heapify(heap, i, 0);
    }

    memcpy(arr, heap, sizeof(int) * n);

    DeAlloc(heap);
}

// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
void merge(int pData[], int l, int m, int r)
{
    int n1 = m - l + 1;
    int n2 = r - m;

    // allocate memory using alloc
    int *L = (int *)Alloc(sizeof(int) * n1);
    int *R = (int *)Alloc(sizeof(int) * n2);

    for (int i = 0; i < n1; i++)
        L[i] = pData[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = pData[m + 1 + j];

    int i = 0;
    int j = 0;
    int k = l; 
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            pData[k] = L[i];
            i++;
        } else {
            pData[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        pData[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        pData[k] = R[j];
        j++;
        k++;
    }

    // deallocate memory using dealloc
    DeAlloc(L);
    DeAlloc(R);
}

// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
void mergeSort(int pData[], int l, int r)
{
    if (l < r) {
        int m = l + (r - l) / 2;
        mergeSort(pData, l, m);
        mergeSort(pData, m + 1, r);
        merge(pData, l, m, r);
    }
}

// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int arr[], int n) {
    int *temp = (int *)Alloc(sizeof(int) * n);
    memcpy(temp, arr, sizeof(int) * n);

    int i, key, j;
    for (i = 1; i < n; i++) {
        key = temp[i];
        j = i - 1;

        while (j >= 0 && temp[j] > key) {
            temp[j + 1] = temp[j];
            j = j - 1;
        }
        temp[j + 1] = key;
    }


    memcpy(arr, temp, sizeof(int) * n);

    DeAlloc(temp);
}

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int arr[], int n) {
    int *temp = (int *)Alloc(sizeof(int) * n);

    memcpy(temp, arr, sizeof(int) * n);

    int i, j;
    for (i = 0; i < n - 1; i++)

        for (j = 0; j < n - i - 1; j++)
            if (temp[j] > temp[j + 1]) {
                int tempVal = temp[j];
                temp[j] = temp[j + 1];
                temp[j + 1] = tempVal;
            }

    memcpy(arr, temp, sizeof(int) * n);

    DeAlloc(temp);
}

// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int arr[], int n) {
    int *temp = (int *)Alloc(sizeof(int) * n);

    memcpy(temp, arr, sizeof(int) * n);

    int i, j, min_idx;

    for (i = 0; i < n - 1; i++) {

        min_idx = i;
        for (j = i + 1; j < n; j++)
            if (temp[j] < temp[min_idx])
                min_idx = j;

        int tempVal = temp[min_idx];
        temp[min_idx] = temp[i];
        temp[i] = tempVal;
    }

    memcpy(arr, temp, sizeof(int) * n);

    DeAlloc(temp);
}

// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
  FILE* inFile = fopen(inputFileName,"r");
  int dataSz = 0;
  *ppData = NULL;

  if (inFile)
  {
    fscanf(inFile,"%d\n",&dataSz);
    *ppData = (int *)Alloc(sizeof(int) * dataSz);
    for (int i = 0; i < dataSz; i++) {
            fscanf(inFile, "%d", &((*ppData)[i]));
        }
        fclose(inFile);
  }

  return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

int main(void)
{
	clock_t start, end;
	int i;
    double cpu_time_used;
	char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};
	
	for (i=0;i<3;++i)
	{
		int *pDataSrc, *pDataCopy;
		int dataSz = parseData(fileNames[i], &pDataSrc);
		
		if (dataSz <= 0)
			continue;
		
		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);
	
		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");
		
		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

                printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}
	
}
