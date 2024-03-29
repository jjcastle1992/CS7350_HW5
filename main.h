/*********************
Name: James Castle
CS 7350 HW5 Graph Project
Purpose: Interface for my driver file that will be testing my graph object.
**********************/

#ifndef HW5_CODE_MAIN_H
#define HW5_CODE_MAIN_H

#include <iostream>
#include <string>
#include <sstream> // for file parsing
#include <bits/stdc++.h>
#include <chrono>
#include <cstdlib> // For function calls in randomRangeGen

struct Node {
    int nodeID;
    Node *nextNode;
    Node *prevNode; // to get the node previous weight
    int nextNodeWeight;
    int numDegrees; // to detect Euler tours at head of each list element
};

class Graph {
public:
    //Constructor
    Graph();

    //Destructor - not going to worry about for now

    //Getters <- breaking encapsulation by not using for attributes, but not trying to be robust.
//    int numNodes();
//    int numEdges();

    // Setters  <- breaking encapsulation by not using for attributes, but not trying to be robust
    void setVerbose(bool); // Set debugging prints on or off
    int getNodeCount(); // Gets nodeCount

    //Vertex & Edge Add/Removal
    void addNode(int); // Takes in a Node ID
    void addEdge(int, int, int); //Takes in starting nodeID, destination Node ID, and weight,
    void addUndirectedEdge(int, int, int); // Calls addEdge twice - one for each direction

    // Printing and utilities
    void displayGraph(); // Displays the entire graph in adjacency list format.
    void clearGraph(); // Deletes Adjacency Matrix
    static bool sortcol(const std::vector<int>&, const std::vector<int>&);

    // Graph Creation Methods
    void createRandomGraph(int, int, bool, unsigned int); // Creates a random graph that is connected, but may undirected or directed.
                                            // Can use random weights or weight of 1.


private:
    int findNode (int); // looks in our adjacency list to find a node. Accepts an int node ID, and a bool for verbose
    std::vector<Node*> adjacencyList; // Our graph in adjacency list format
    int nodeCount; // counts number of nodes in graph (should be == to size of vector)
    int edgeCount; // Counts each unique edge between nodes (undirected edges get counted once)
    int maxDegree; // Counts the maximum degree of any vertex in the graph
    std::vector<int> oddDegreeNodes; // List containing all nodes of odd degree (Determines Euler Tour Truth)
    bool verbose; // turns on/off debug prints
    int maxWeight; // max edge weight in the graph
};

int randomRangeGen(int, int, unsigned int);

#endif //HW5_CODE_MAIN_H
