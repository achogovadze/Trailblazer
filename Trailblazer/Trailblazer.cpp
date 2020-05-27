/******************************************************************************
 * File: Trailblazer.cpp
 *
 * Implementation of the graph algorithms that comprise the Trailblazer
 * assignment.
 */

#include "Trailblazer.h"
#include "TrailblazerGraphics.h"
#include "TrailblazerTypes.h"
#include "TrailblazerPQueue.h"
#include "random.h"
using namespace std;

/* Function: shortestPath
 * 
 * Finds the shortest path between the locations given by start and end in the
 * specified world.	 The cost of moving from one edge to the next is specified
 * by the given cost function.	The resulting path is then returned as a
 * Vector<Loc> containing the locations to visit in the order in which they
 * would be visited.	If no path is found, this function should report an
 * error.
 *
 * In Part Two of this assignment, you will need to add an additional parameter
 * to this function that represents the heuristic to use while performing the
 * search.  Make sure to update both this implementation prototype and the
 * function prototype in Trailblazer.h.
 */

// this method finds all neighbours of the given location
Set<Loc> visitNeighbours (Grid<double>& world, Loc cur){
	Set<Loc> result;
	int startLocCol = cur.col - 1;
	int startLocRow = cur.row - 1;
	for (int i = startLocRow; i <= startLocRow + 2; i++){
		for (int k = startLocCol; k <= startLocCol + 2; k++){
			if (i == cur.row && k == cur.col) continue;
			if (world.inBounds(i,k)){
				Loc loc = makeLoc(i, k);
				result += loc;
			}
		}
	}
	return result;
}

Vector<Loc>
shortestPath(Loc start,
             Loc end,
             Grid<double>& world,
             double costFn(Loc from, Loc to, Grid<double>& world),
			 double heuristic(Loc start, Loc end, Grid<double>& world)
			 ) {
	// creating structures
	Vector<Loc> result;
	TrailblazerPQueue<Loc> queue;
	Map<Loc, Loc> parents;
	Map<Loc, Color> colors;
	Map<Loc, double> dist;

	// initializing first Loc
	Loc finish = end;
	colors[start] = YELLOW;
	colorCell(world, start, YELLOW);
	dist[start] = 0;
	queue.enqueue(start, heuristic(start, end, world));

	// searching shortest path
	while (!queue.isEmpty()){
		Loc cur = queue.dequeueMin();
		colors[cur] = GREEN;
		colorCell(world, cur, GREEN);
		double curDist = dist[cur];
		if (cur == finish){
			finish = cur;
			break;
		}
		Set<Loc> neighbours = visitNeighbours(world, cur);
		foreach (Loc n in neighbours){
			double L = costFn(cur, n, world);
			if (!colors.containsKey(n)){
				colors[n] = YELLOW;
				colorCell(world, n, YELLOW);
				dist[n] = curDist + L;
				parents[n] = cur;
				queue.enqueue(n, dist[n] + heuristic(n, end, world));
			} else if(colors[n] == YELLOW && dist[n] > curDist + L){
				dist[n] = curDist + L;
				parents[n] = cur;
				queue.decreaseKey(n, dist[n] + heuristic(n, end, world));
			}
		}
	}

	// fills final result vector with locations
	while (finish != start){
		result.insert(0, finish);
		finish = parents[finish];
	}
	result.insert(0, start);
	return result;
}

// this method creates graph in a form of a grid
void setupGrid(Grid<Loc>& grid, int numRows, int numCols){
	grid.resize(numRows, numCols);
	for (int i = 0; i < numRows; i++){
		for (int k = 0; k < numCols; k++){
			grid[i][k] = makeLoc(i, k);
		}
	}
}

// this method adds all edges of the graph to a priority queue
void enqueueEdges(Grid<Loc>& grid, TrailblazerPQueue<Edge>& queue){
	for (int i = 0; i < grid.nRows; i++){
		for (int c = 0; c < grid.nCols; c++){
			if (grid.inBounds(i, c + 1)){
				Loc start = grid[i][c];
				Loc end = grid[i][c + 1];
				Edge newEdge = makeEdge(start, end);
				queue.enqueue(newEdge, randomReal(0, 1));
			}
			if (grid.inBounds(i + 1, c)){
				Loc start = grid[i][c];
				Loc end = grid[i + 1][c];
				Edge newEdge = makeEdge(start, end);
				queue.enqueue(newEdge, randomReal(0, 1));
			}
		}
	}
}

// this method fills map with key of a location of node and with value a set of nodes connected to the original node
void putClusters(Grid<Loc>& grid, Map<Loc, Set<Loc> >& clusters){
	for (int i = 0; i < grid.nRows; i++){
		for (int c = 0; c < grid.nCols; c++){
			Loc cur = grid[i][c];
			Set<Loc> locations;
			locations.add(cur);
			clusters.put(cur, locations);
		}
	}
}

// this method merges two clusters
void merge (Map<Loc, Set<Loc> >& clusters, Edge e){
	Loc toChange;
	Loc toDelete;
	foreach (Loc l in clusters){
		if (clusters[l].contains(e.start))
			toChange = l;
		if (clusters[l].contains(e.end))
			toDelete = l;
	}
	foreach (Loc l in clusters[toDelete]){
		clusters[toChange] += l;
	}
	clusters.remove(toDelete);
}

// this method checkes if two locations are in a same cluster
bool sameCluster(Loc start, Loc end, Map<Loc, Set<Loc> >& clusters){
	foreach(Loc l in clusters){
		if (clusters[l].contains(start) && clusters[l].contains(end))
			return true;
	}
	return false;
}

Set<Edge> createMaze(int numRows, int numCols) {
	Set<Edge> result;
	Grid<Loc> grid;
	TrailblazerPQueue<Edge> queue;
	Map<Loc, Set<Loc> > clusters;

	setupGrid(grid, numRows, numCols);
	enqueueEdges(grid, queue);
	putClusters(grid, clusters);

	while (clusters.size() >= 2){
		Edge cur = queue.dequeueMin();
		if (!sameCluster(cur.start, cur.end, clusters)){
			merge(clusters, cur);
			result += cur;
		}
	}
	return result;
}
