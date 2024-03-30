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
    maxWeight = 0;
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

                if(edgeWeight > this->maxWeight) {
                    maxWeight = edgeWeight;
                }

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
              this->maxDegree << "\nMax Weight: " << this->maxWeight << std::endl;
    std::cout << "Graph adjacency list printed" << std::endl;
}

void Graph::clearGraph() {
    this->adjacencyList.clear();
    this->nodeCount = 0;
    this->edgeCount = 0;
}

void Graph::writeToFile(const std::string& filename) {
    std::vector<Node*> graphToWrite = this->adjacencyList;
    std::ofstream outfile;
//    std::cout << filename;
    outfile.open(filename);

    if(outfile.fail()){ // File could not be opened
        std::cerr << "ERROR: File could not be opened" << std::endl;
        std::exit(1);
    }
    outfile << "Start Node, Destination Node, Edge Weight" << std::endl;

    Node *currentNode;
    for(int i = 0; i < this->adjacencyList.size(); i++){
        // Check to see if node has an edge and write out
        int startId = -1;
        currentNode = this->adjacencyList[i];
        startId = currentNode->nodeID;

        while(currentNode->nextNode != nullptr){
            int destId = -1;
            int weight = -1;
            destId = currentNode->nextNode->nodeID;
            weight = currentNode->nextNodeWeight;
            // If all three are valid (not -1), then write to file
            if ((startId != -1) && (destId != -1) && (weight != -1)){
                outfile << startId << "," << destId << "," << weight << std::endl;
            }
            currentNode = currentNode->nextNode;
        }
    }
    outfile.close();
}

void Graph::createRandomGraph(int numVerts, int numEdges, bool randomWeights, unsigned int seed = 42, double skew = 1) {
    /* Accepts a number of verts and a number of edges and randomly creates edges between 2 random vertices (that may or
    may not be undirected). The graph may or may not be connected, and can have anywhere from 0 nodes with 0 edges to V
    vertices with (v * (v - 1))/2 edges <- I.e. a complete graph. There is an option to use Random weights or a weight
     of 1 for all edges.
    */
    bool noSetSeed = false;
    if (seed == 42){
        noSetSeed = true; // no seed given, go random
    }
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
            if(noSetSeed){
                startNodeId = skewRandomRangeGen((numVerts - 1), 0, 42, skew);
                endNodeId = skewRandomRangeGen((numVerts - 1), 0, 42, skew);

                // Determine if using random weights or not
                if(randomWeights){
                    weight = randomRangeGen(100, 1, seed);  // Only edge choices should be skewed
                }
                else{
                    weight = 1;
                }
            }
            else {
                startNodeId = skewRandomRangeGen((numVerts - 1), 0, seed, skew);
                endNodeId = skewRandomRangeGen((numVerts - 1), 0, seed + 1, skew);
                seed++;
                // Determine if using random weights or not
                if(randomWeights){
                    weight = randomRangeGen(100, 1, seed);  // Only edge choices should be skewed
                }
                else{
                    weight = 1;
                }
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
    int random;

    // Random pathway
    if (seed == 42) {
        random = startRange + (rand() % ((endRange - startRange) + 1));
    }
        // Set seed pathway
    else {
        // Modified with ChatGPT to take in a seed.
        // Initialize the random number generator with the provided seed
        std::mt19937 gen(seed);
        // Retrieve a random number between startRange and EndRange
        random = startRange + (gen() % ((endRange - startRange) + 1));
    }
    return random;
}

// Function to generate skewed random numbers
int skewRandomRangeGen(int endRange, int startRange = 0, unsigned int seed = 42, double skew = 1.0) {
    // From ChatGPT
    // Initialize the random number generator with the provided seed or a random seed
    std::mt19937 gen;
    if (seed == 42) {
        // Generate a random seed based on the current time
        seed = static_cast<unsigned int>(std::time(nullptr));
        gen.seed(seed);
    } else {
        // Initialize the random number generator with the provided seed
        gen.seed(seed);
    }

    // Generate a random number between 0 and 1
    std::uniform_real_distribution<double> dist(0, 1);
    double randNum = dist(gen);

    // Apply skewness to the random number
    double skewedRand = pow(randNum, skew);

    // Scale the skewed random number to fit the range
    int random = startRange + static_cast<int>((endRange - startRange + 1) * skewedRand);

    return random;
}

int main() {
    unsigned int seed = 4002; // NOTE: to set the seed manually, choose any value OTHER THAN 42.

    std::srand(time(nullptr));

    int maxEdges = 0;
    int stepSize = 10;
    int maxSize = 100;
    std::vector<std::vector<long long>>  randomWeightTiming;  // [NumNodes][Q1 Time (Graph read)]]
    std::vector<std::vector<long long>>  weightOneTiming;  // [NumNodes][Q1 Time (Graph read)]]

    // Test with random weights
    std::cout << "*************Graphs with Random Weights*************" << std::endl;
    for (int i = 10; i <= maxSize; i += stepSize){
        Graph randomWeights;
        maxEdges = ((i * (i - 1)) / 2);
        int edges = skewRandomRangeGen(maxEdges, 0, seed + 1);

        // Q1 Timing
        auto q1start = std::chrono::high_resolution_clock::now();  // start timer
        randomWeights.createRandomGraph(i, edges, true, seed, 1);  // Note a seed of 42 will be randomized
        randomWeights.displayGraph();
        auto q1stop = std::chrono::high_resolution_clock::now(); // stop timer
        auto q1duration = std::chrono::duration_cast<std::chrono::microseconds>(q1stop - q1start);
        randomWeightTiming.push_back({randomWeights.getNodeCount(), q1duration.count()});

        //Write to file
        std::string filename = "randomGraphRandomWeight_" + std::to_string(i) + ".csv";
        randomWeights.writeToFile(filename);

        randomWeights.clearGraph();
        std::cout << "\nGraph Size: " << i << std::endl;
    }

     // Test without random weights
    std::cout << "*************Graphs with Weight 1*************" << std::endl;
    for (int i = 10; i <= maxSize; i += stepSize){
        Graph weight1;
        maxEdges = ((i * (i - 1)) / 2);
        int edges = randomRangeGen(maxEdges, 0, seed + i);
        // Q1 Timing
        auto q1start = std::chrono::high_resolution_clock::now();  // start timer
        weight1.createRandomGraph(i, edges, false, seed, 1);  // Note a seed of 42 will be randomized
        auto q1stop = std::chrono::high_resolution_clock::now(); // stop timer
        auto q1duration = std::chrono::duration_cast<std::chrono::microseconds>(q1stop - q1start);
        weightOneTiming.push_back({weight1.getNodeCount(), q1duration.count()});

        // Write to file
        std::string filename = "randomWeight1_" + std::to_string(i) + ".csv";
        weight1.writeToFile(filename);

        weight1.clearGraph();
        std::cout << "\nGraph Size: " << i << std::endl;
    }

    // Display Timing
    std::cout << "\nRandom Weight Timing Stats:";
    for (int i = 0; i < randomWeightTiming.size(); i++) {
      std::cout << "\n\tNodeCount: " << randomWeightTiming[i][0] <<
      "\n\tQ1 (Create) Timing: " << randomWeightTiming[i][1] << " Microseconds" << std::endl;
    }

    std::cout << "\nWeight 1 Timing Stats:";
    for (int i = 0; i < weightOneTiming.size(); i++) {
        std::cout << "\n\tNodeCount: " << weightOneTiming[i][0] <<
                  "\n\tQ1 (Create) Timing: " << weightOneTiming[i][1] << " Microseconds" << std::endl;
    }

    std::cout << "\ndone" << std::endl;
    return 0;
}