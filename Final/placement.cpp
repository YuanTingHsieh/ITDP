#include "evaluate.h"
#include <iostream>
#include <vector>
#include <list>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "RBTree.h"
using namespace std;

void IntDest(void* a) {
  free((int*)a);
}

int IntComp(const void* a,const void* b) {
  if(*(int*)a > *(int*)b) return(1);
  if(*(int*)a < *(int*)b) return(-1);
  return(0);
}

void IntPrint(const void* a) {
  printf("%i",*(int*)a);
}

void InfoPrint(void* a) {
	printf("%i",*(int*)a);
}

void InfoDest(void *a){
  ;
}

void sizeInfoPrint(void* a) {
	rb_red_blk_tree* distTree = (rb_red_blk_tree*)a;
	printf("\n----------------------------------------------------------------------\n");
	RBTreePrint(distTree);
}

bool circuit::isOverlapped(int cellIndex) {
	cell* theCell = &cells[cellIndex];
	int rowIndex = (theCell->y_coord - rows[0].origY) / (int)rowHeight;
	int siteOffset = (theCell->x_coord - rows[0].origX) / rows[0].stepX;
	int siteIndex = rowIndex * siteMappingWidth + siteOffset;
	int cellSiteNum = theCell->width / rows[0].stepX;
	bool result = false;
	for(int i = 0; i < cellSiteNum && i < siteMappingWidth; i++) {
		if(siteMapping[rowIndex][siteOffset + i] > 1)
			result = true;
	}
	return result;
}

/*void circuit::removeOverlapping() {
	for(unsigned i = 0; i < cells.size(); i++) {
		cell* theCell = &cells[i];
		if(!theCell->isFixed && isOverlapped(i)) {
			int rowIndex = (theCell->y_coord - rows[0].origY) / (int)rowHeight;
			int siteOffset = (theCell->x_coord - rows[0].origX) / rows[0].stepX;
			int siteIndex = rowIndex * siteMappingWidth + siteOffset;
			int cellSiteNum = theCell->width / rows[0].stepX;
			
			for(int j = 0; j < cellSiteNum; j++)
				siteMapping[rowIndex][siteOffset + j]--;
			bool find = false;
			for(int offsetY = 0; !find && offsetY < 3 && rowIndex + offsetY < rows.size(); offsetY++) {
				for(int offsetX = 0; !find && offsetX + cellSiteNum <= siteMappingWidth; offsetX++) {
					for(int j = 0; j < cellSiteNum; j++) {
						if(siteMapping[rowIndex + offsetY][offsetX + j] > 0) {
							find = true;
							break;
						}
					}
					for(int j = 0; j < cellSiteNum; j++)
						siteMapping[rowIndex + offsetY][offsetX + j]++;
					theCell->y_coord = (rowIndex + offsetY) * (int)rowHeight + rows[0].origY;
					theCell->x_coord = offsetX * rows[0].stepX + rows[0].origX;
				}
			}
		}
	}
}*/

void circuit::removeOverlapping() {
	printf("---------------------removeOverlapping()--------------------------------\n");
	for(unsigned i = 0; i < cells.size(); i++) {
		cell* theCell = &cells[i];
		if(!theCell->isFixed && isOverlapped(i)) { // move the cell to remove overlapping
			int rowIndex = (theCell->y_coord - rows[0].origY) / (int)rowHeight;
			int siteOffset = (theCell->x_coord - rows[0].origX) / rows[0].stepX;
			int siteIndex = rowIndex * siteMappingWidth + siteOffset;
			int cellSiteNum = theCell->width / rows[0].stepX;
			// choose which distTree to insert
			rb_red_blk_node* sizeNode = findDistTree(sizeTree, &cellSiteNum, 1);
			rb_red_blk_tree* distTree = (rb_red_blk_tree*)sizeNode->info;
			rb_red_blk_node* freeSiteNode = findFreeSites(distTree, &siteIndex);
			while(freeSiteNode == NULL) {
				sizeNode = TreeSuccessor(sizeTree, sizeNode);
				distTree = (rb_red_blk_tree*)sizeNode->info;
				freeSiteNode = findFreeSites(distTree, &siteIndex); 
			}
			moveCell(i, *(int*)freeSiteNode->key, *(int*)freeSiteNode->info);
			free(freeSiteNode);
		}
	}
	for(int i = 0; i < rows.size(); i++) {
		for(int j = 0; j < siteMappingWidth; j++) {
			printf("%d ", siteMapping[i][j]);    
		}
		printf("\n");
	}
}


void circuit::buildDistTreesQQQ() {
	for(int i = 0; i < rows.size(); i++) {
		int count = 0;
		for(int j = 0; j < siteMappingWidth; j++) {
			
			if(siteMapping[i][j] == 0)
				count++;
			else if(count != 0) { // insert distTree
				// choose which distTree to insert
				rb_red_blk_node* sizeNode = findDistTree(sizeTree, &count, 0);
				rb_red_blk_tree* distTree = (rb_red_blk_tree*)(sizeNode->info);
				int* siteIndex = (int*)malloc(sizeof(int));
				int* numSites = (int*)malloc(sizeof(int));
				*siteIndex = i * siteMappingWidth + j - count;
				*numSites = count;
				RBTreeInsert(distTree, siteIndex, numSites);
				count = 0;
			}
		}
		if(count != 0) {
			// choose which distTree to insert
			rb_red_blk_node* sizeNode = findDistTree(sizeTree, &count, 0);
			rb_red_blk_tree* distTree = (rb_red_blk_tree*)(sizeNode->info);
			int* siteIndex = (int*)malloc(sizeof(int));
			int* numSites = (int*)malloc(sizeof(int));
			*siteIndex = (i + 1) * siteMappingWidth - count;
			*numSites = count;
			RBTreeInsert(distTree, siteIndex, numSites);
		}
	}
}

void circuit::buildSizeTree() {
	sizeTree = RBTreeCreate(IntComp, IntDest, InfoDest, IntPrint, sizeInfoPrint);
	// default size 1 distTree
	int* initSize = (int*)malloc(sizeof(int));
	*initSize = 1;
	rb_red_blk_node* node = RBTreeInsert(sizeTree, initSize, 0);
	node->info = RBTreeCreate(IntComp, IntDest, InfoDest, IntPrint, InfoPrint);
	
	for(unsigned i = 0; i < cells.size(); i++) {
		int cellSiteNum = cells[i].width / rows[0].stepX;
		int* size = (int*)malloc(sizeof(int));
		*size = cellSiteNum;
		if(!RBExactQuery(sizeTree, size)) {
			node = RBTreeInsert(sizeTree, size, 0);
			node->info =  RBTreeCreate(IntComp, IntDest, InfoDest, IntPrint, InfoPrint);
		}
		else
			free(size);
	}
	
	//RBTreePrint(sizeTree);
}

void circuit::buildSiteMapping() {
	siteMapping = new int*[rows.size()];
	for(int i = 0; i < rows.size(); i++) {
		siteMapping[i] = new int[siteMappingWidth];
		for(int j = 0; j < siteMappingWidth; j++)
			siteMapping[i][j] = 0;
		for(int j = rows[i].numSites; j < siteMappingWidth; j++)
			siteMapping[i][j] = 1;
	}
	for(vector<cell>::iterator it = cells.begin(); it != cells.end(); it++) {
		
		int rowIndex = (it->y_coord - rows[0].origY) / (int)rowHeight;
		int siteOffset = (it->x_coord - rows[0].origX) / rows[0].stepX;
		int cellSiteNum = it->width / rows[0].stepX;
		//printf("cell %s = (%d, %d, %d)\n", it->name.c_str(), rowIndex, siteOffset, cellSiteNum);
		for(int i = 0; i < cellSiteNum && siteOffset + i < siteMappingWidth; i++)
				siteMapping[rowIndex][siteOffset + i]++;
	}
	printf("----------------------buildSiteMapping()----------------------\n");
	for(int i = 0; i < rows.size(); i++) {
		for(int j = 0; j < siteMappingWidth; j++) {
			printf("%d ", siteMapping[i][j]);    
		}
		printf("\n");
	}
} 
void circuit:: moveCell(int cellIndex, int freeSiteIndex, int freeSiteSize) {
	cell* theCell = &cells[cellIndex];
	int rowIndex = (theCell->y_coord - rows[0].origY) / (int)rowHeight;
	int siteOffset = (theCell->x_coord - rows[0].origX) / rows[0].stepX;
	int siteIndex = rowIndex * siteMappingWidth + siteOffset;
	int cellSiteNum = theCell->width / rows[0].stepX;
	
	int freeSiteRowIndex = freeSiteIndex / siteMappingWidth;
	int freeSiteOffset = freeSiteIndex % siteMappingWidth;
	theCell->y_coord = freeSiteRowIndex * (int)rowHeight + rows[0].origY;
	theCell->x_coord = freeSiteOffset * rows[0].stepX + rows[0].origX;
	
	for(int i = 0; i < cellSiteNum && siteOffset + i < siteMappingWidth; i++) {
		siteMapping[rowIndex][siteOffset + i]--;
		siteMapping[freeSiteRowIndex][freeSiteOffset + i]++;
	}
}

void circuit::traverseCells() {
	for(list<unsigned>::iterator it = traverseList.begin(); it != traverseList.end(); it++) {
		unsigned cellIndex = *it;
		//moveCell(cellIndex);
	}
}
void circuit::buildTraverseOrder() {
	
	list<unsigned> pinsQueue;
	// initial
	for(vector<cell>::iterator it = cells.begin(); it != cells.end(); it++) {
		it->marked = false;
	}
	
	// traverse
	for(unsigned int i = 0; i < PIs.size(); i++) {
		unsigned pinIndex = PIs[i];
		pinsQueue.push_back(pinIndex);
	}
	while(!pinsQueue.empty()) {
		unsigned pinIndex = pinsQueue.front();
		pinsQueue.pop_front();
		pin* thePin = &pins[pinIndex];
		//cout << "port: " << thePin->name << " is visited." << endl;
		if(thePin->portType == 2 || thePin->portType == 1) {// primary Input or output port
			unsigned netIndex = thePin->net;
			net* theNet = &nets[netIndex];
			for(unsigned int i = 0; i < theNet->sinks.size(); i++) {
				unsigned sinkPinIndex = theNet->sinks[i];
				pinsQueue.push_back(sinkPinIndex);
			}
		}
		else if(thePin->portType == 0) { // input port
			// add the cell into list & push output pins
			unsigned cellIndex = thePin->owner;
			cell* theCell = &cells[cellIndex];
			if(theCell->marked)
				continue;
			// visit the cell
			traverseList.push_front(cellIndex);
			theCell->marked = true;
			// push output pins
			for(map<string, unsigned>::iterator it = theCell->ports.begin(); it != theCell->ports.end(); it++) {
				pinIndex = it->second;
				thePin = &pins[pinIndex];
				if(thePin->portType == 1) // output port
					pinsQueue.push_back(pinIndex);
			}
		}
	}
}

void circuit::outputDEF(ofstream &os) {
	os << "VERSION " << DEFVersion << " ;" << endl;
	os << "DESIGN " << design_name << " ;" << endl << endl;
	os << "COMPONENTS " << cells.size() << " ;" << endl;
	for(unsigned i = 0; i < cells.size(); i++) {
		cell* theCell = &cells[i];
		macro* theMacro = &macros[theCell->type];
		os << "   - " << theCell->name << " " << theMacro->name << endl;
		if(theCell->isFixed)
			os << "      + " << "FIXED ( " << theCell->x_coord << " " << theCell->y_coord << " ) " << theCell->cellorient << " N " << endl;
		else
			os << "      + " << "PLACED ( " << theCell->x_coord << " " << theCell->y_coord << " ) " << theCell->cellorient << " N " << endl;
	}
	os << "END COMPONENTS" << endl << endl;
	os << "END DESIGN" << endl;
}






	
