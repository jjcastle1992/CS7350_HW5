/*********************
Name: James Castle
CS 7350 HW5 Graph Project
Purpose: Implementation for my driver file that will be testing my graph object.
**********************/

#include "main.h"

// Concepts borrowed from
// https://www.programiz.com/dsa/graph-adjacency-list
// https://www.softwaretestinghelp.com/graph-implementation-cpp/

// Constructor
Graph::Graph() {
    nodeCount = 0;
    edgeCount = 0;
    maxDegree = 0;
    verbose = false;
    this->oddDegreeNodes;
}

void Graph::setVerbose(bool beVerbose) {
    this->verbose = beVerbose;
}

int Graph::getNodeCount() {
    return this->nodeCount;
}

void Graph::addNode(int id) {
    Node *newNode = new Node;
    newNode->nodeID = id;
    newNode->nextNode = nullptr;
    newNode->prevNode = nullptr;
    newNode->nextNodeWeight = -1; // Signals there are no connected nodes
    newNode->numDegrees = 0;
    if(this->adjacencyList.empty()){
        this->adjacencyList.push_back(newNode);
    }
    if(this->adjacencyList.size() < (id + 1)) {
        this->adjacencyList.resize(id + 1);
    }
//    adjacencyList.push_back(newNode);
    this->adjacencyList[id] = newNode;
    this->nodeCount++;
}

void Graph::addEdge(int startId, int destId, int edgeWeight = 1) {
    //Check to see if the start and dest nodes are not the same
    if(startId != destId){
        // validate that start and destination nodes exist
        startId = findNode(startId);
        destId = findNode(destId);
        bool preExistingEdge = false;
        int localMaxDegree = 0;  //Prep to check if max degree is greater

        if((startId != -1) && (destId != -1)) {
            // Validate edge does not already exist
            //only check start LL for connection to dest
            Node *headNode = adjacencyList[startId];
            Node *currentNode = adjacencyList[startId];
            Node *temp = adjacencyList[destId];
            Node *newEdge = new Node;
            newEdge->nodeID = temp->nodeID;
            newEdge->nextNode = nullptr;
            newEdge->prevNode = nullptr;
            newEdge->nextNodeWeight = -1;

            while(currentNode->nextNode != nullptr) {
                if(currentNode->nextNode->nodeID == destId){
                    if(this->verbose){
                        std::cout << "Edge Already exists" << startId << " -> " << destId << std::endl;
                    }
                    preExistingEdge = true;
                }
                currentNode = currentNode->nextNode;
                localMaxDegree++;
            }
            if (!preExistingEdge) { // if the edge does not exist, add it.
                currentNode->nextNode = newEdge;
                newEdge->nextNode = nullptr;
                newEdge->prevNode = currentNode;
                currentNode->nextNodeWeight = edgeWeight;
                localMaxDegree++;
                headNode->numDegrees += 1;
                if(headNode->numDegrees % 2 != 0){ // Checks to see if the node has odd degrees
                    //Odd number of degrees, push on oddDegreeNodes list
                    this->oddDegreeNodes.push_back(headNode->nodeID);
                }
                else {
                    // Find and remove ID from oddDegreeNodes list;
                    // Used this source: https://www.shiksha.com/online-courses/articles
                    // /erasing-elements-from-a-vector-in-c/
                    for(int i = 0; i < this->oddDegreeNodes.size(); i++){
                        if(this->oddDegreeNodes[i] == headNode->nodeID){
                            this->oddDegreeNodes.erase(this->oddDegreeNodes.begin() + i);
                            break;
                        }
                    }
                }

                // check for previous directed edge (as undirected edge should be counted once)
                Node *undirCheck = this->adjacencyList[destId];
                bool undirected = false;
                if (undirCheck != nullptr){
                    while(undirCheck->nextNode != nullptr){
                        if(undirCheck->nextNode->nodeID == startId) {
                            undirected = true;
                        }
                        undirCheck = undirCheck->nextNode;
                    }
                }
                if(!undirected) {
                    this->edgeCount++;
                }
            }
        }
        else{
            std::cout << "ERROR: One or more passed in edges are invalid" << std::endl;
        }
        if(localMaxDegree > this->maxDegree) {
            this->maxDegree = localMaxDegree;
        }
    }
}

void Graph::addUndirectedEdge(int startId, int destId, int weight) {
    addEdge(startId, destId, weight);
    addEdge(destId, startId, weight);
}

int Graph::findNode(int targetNodeId) {
    int nodeId = -1; // indicates not found
    //Determine if graph is empty
    if(adjacencyList.empty()){
        if(verbose){
            std::cout << "Graph Empty" << std::endl;
        }
    }
    else{
        for (int i = 0; i < adjacencyList.size(); i++) {
            if(adjacencyList[i] == nullptr){
                continue;
            }
            int currentNodeId = adjacencyList[i]->nodeID;
            if (targetNodeId == currentNodeId) {
                nodeId = currentNodeId;
            }
        }
    }
    if(verbose){
        if(nodeId == -1) {
            std::cout << "node not found" << std::endl;
        }
        else {
            std::cout << "node found" << std::endl;
        }
    }

    return nodeId;
}
void Graph::displayGraph() {
    std::vector<Node*> graphToPrint = this->adjacencyList;
    std::cout << "\nList of nodes in our adjacency list" << std::endl;
    for (int i = 0; i < graphToPrint.size(); i++ ) {
        int degrees = 0;
        bool firstPrint = true;
        // Print node at Root
        std::cout << "Node "<< graphToPrint[i]->nodeID << ": ";
        Node *currentNode = graphToPrint[i];
        while(currentNode->nextNode != nullptr){
            if(firstPrint){
                std::cout << graphToPrint[i]->nodeID;
            }
            currentNode = currentNode->nextNode;
            std::cout << " --> " << currentNode->nodeID << "(" << currentNode->prevNode->nextNodeWeight << ")";
            degrees++;
            firstPrint = false;
        }
        if(degrees == 0) {
            std::cout << "EMPTY";
        }
        std::cout << std::endl;
    }
    std::cout << "# nodes: " << this->nodeCount << "\n# edges: " << this->edgeCount << "\nMax Degrees: " <<
              this->maxDegree << "\nEuler Tour: " << (this->oddDegreeNodes.empty() ? "true":"false") << std::endl;
    std::cout << "Graph adjacency list printed" << std::endl;
}

void Graph::clearGraph() {
    this->adjacencyList.clear();
    this->nodeCount = 0;
    this->edgeCount = 0;
}

void Graph::createRandomGraph(int numVerts, int numEdges, bool randomWeights, unsigned int seed = 42) {
    /* Accepts a number of verts and a number of edges and randomly creates edges between 2 random vertices (that may or
    may not be undirected). The graph may or may not be connected, and can have anywhere from 0 nodes with 0 edges to V
    vertices with (v * (v - 1))/2 edges <- I.e. a complete graph. There is an option to use Random weights or a weight
     of 1 for all edges.
    */


    // Create adjacency list
    for (int i = 0; i < numVerts; i++) {
        this->addNode(i);
    }

    // Edge Creation and validation
    int maxEdges = ((numVerts * (numVerts - 1)) / 2);
    if(numEdges > maxEdges){
        std::cout << "Number of requested edges (" << numEdges << " ) exceeds total maximum possible edges ("
                  << maxEdges << "). Replacing edge request total with MaxEdges of: " << maxEdges << std::endl;
        numEdges = maxEdges;
    }

    while (this->edgeCount < numEdges){
        // Use Random to pick a start vert and dest vert between 0 and V - 1 where V - 1 is the last element in the
        // adj list
        int startNodeId = -1;
        int endNodeId = -1;
        int weight = -1;

        // Pick our random nodes to make an edge with
        do {

            startNodeId = randomRangeGen((numVerts - 1), 0, seed);
            endNodeId = randomRangeGen((numVerts - 1), 0, seed + 1);
            seed++;
            // Determine if using random weights or not
            if(randomWeights){
                weight = randomRangeGen(100, 1, seed);
            }
            else{
                weight = 1;
            }
        }
        while(startNodeId == endNodeId);  // Ensure nodes are different

        this->addEdge(startNodeId, endNodeId, weight);  // Add edge to adjacency list
    }
}

bool Graph::sortcol(const std::vector<int> &v1, const std::vector<int> &v2) {
    // From https://www.geeksforgeeks.org/sorting-2d-vector-in-c-set-2-in-descending-order-by-row-and-column/
    return v1[2] < v2[2]; // sort by weight (smallest to largest)
}


int randomRangeGen(int endRange, int startRange = 0, unsigned int seed = 42) {
    // General implementation borrowed from:
    // https://www.digitalocean.com/community/tutorials/random-number-generator-c-plus-plus

    // Modified with ChatGPT to take in a seed.
    // Initialize the random number generator with the provided seed
    std::mt19937 gen(seed);
    // Retrieve a random number between startRange and EndRange
    int random = startRange + (gen() % ((endRange - startRange) + 1));

    return random;
}

// Function to generate skewed random numbers


int main() {

    Graph test;
    test.createRandomGraph(4, 3, true, 45);
    test.displayGraph();
    std::cout << "\ndone" << std::endl;
    return 0;
}