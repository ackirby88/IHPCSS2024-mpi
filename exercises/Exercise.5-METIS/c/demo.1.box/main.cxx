#include <iostream>
#include <fstream>
#include <metis.h>

using namespace std;

// Struct for storing node coordinates
struct Node {
    double x, y;
};

// Struct for storing an element (e.g., triangle) defined by its nodes
struct Element {
    int nodes[3]; // For triangles
};

// Function to write mesh in VTK format
void writeMeshVTK(const char* filename, int numNodes, Node* nodes, int numElems, Element* elems, idx_t* part) {
    ofstream file(filename);

    // Write VTK header
    file << "# vtk DataFile Version 3.0" << endl;
    file << "Partitioned Mesh Example" << endl;
    file << "ASCII" << endl;
    file << "DATASET UNSTRUCTURED_GRID" << endl;

    // Write nodes
    file << "POINTS " << numNodes << " double" << endl;
    for (int i = 0; i < numNodes; ++i) {
        file << nodes[i].x << " " << nodes[i].y << " 0.0" << endl; // Assuming 2D mesh (z=0.0)
    }

    // Write elements
    int totalNodes = 0;
    for (int i = 0; i < numElems; ++i) {
        totalNodes += 1 + elems[i].nodes[2]; // Number of nodes in the current element
    }
    file << "CELLS " << numElems << " " << totalNodes << endl;
    for (int i = 0; i < numElems; ++i) {
        file << "3 "; // Three nodes per triangle (assuming triangles)
        file << elems[i].nodes[0] << " " << elems[i].nodes[1] << " " << elems[i].nodes[2] << endl;
    }

    // Write element types
    file << "CELL_TYPES " << numElems << endl;
    for (int i = 0; i < numElems; ++i) {
        file << "5 "; // VTK_TRIANGLE type
    }
    file << endl;

    // Write partition data as a scalar field
    file << "CELL_DATA " << numElems << endl;
    file << "SCALARS Partition_ID int 1" << endl;
    file << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < numElems; ++i) {
        file << part[i] << endl;
    }

    file.close();
}

int main() {
    // Example mesh (2D unstructured)
    idx_t numNodes = 8;
    Node nodes[numNodes] = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0},
        {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {3.0, 1.0}
    };

    idx_t numElems = 6; // Assuming a simple grid of two triangles
    Element elems[numElems] = {
        {{0, 1, 5}}, {{0, 5, 4}}, {{1, 2, 6}},
        {{1, 6, 5}}, {{2, 3, 7}}, {{2, 7, 6}}
    };

    const int nPi = 3; 
    idx_t *eptr = new idx_t[numElems + 1];
    eptr[0] = 0;
    for(int i = 0 ; i < numElems; i++) eptr[i+1] = eptr[i] + nPi;
    
    idx_t *eind = new idx_t[eptr[numElems]];
    int cnt = 0;
    for (int i = 0 ; i < numElems; i++){
        for (int j = 0; j < nPi; j++){
            eind[cnt] = elems[i].nodes[j];
            cnt++;
        }
    }

    // Convert the mesh to METIS-compatible format
    idx_t ncommon = 2;  // Number of common nodes that two elements must have in order to an edge between them in the dual graph.
    idx_t numParts = 3; // Number of partitions (3-way partitioning)
    
    idx_t *epart = new idx_t[numElems]; // Stores the partition vector of the elements of the mesh
    idx_t *npart = new idx_t[numNodes]; // Stores the partition vector of the nodes of the mesh

    // METIS options
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // Partition the mesh using METIS
    idx_t objval; // METIS output (total edge-cut or communication volume)
    METIS_PartMeshDual(&numElems,   // number of mesh elements
                       &numNodes,// number of mesh nodes
                       eptr,        // array storing mesh node counts per element (scan)
                       eind,        // array storing mesh info
                       NULL,        // vwgt: An array of size ne specifying the weights of the elements.
                       NULL,        // vsize
                       &ncommon,    // number of common nodes two elements need to share a face
                       &numParts,   // number of partitions to generate
                       NULL,        // tpwghts: desired weight for each partition
                       options,     //
                       &objval,     //
                       epart,       // partition number of each element
                       npart);      // partition number of each node
    
    // Output partitioned mesh in VTK format
    writeMeshVTK("partitioned_mesh.vtk", numNodes, nodes, numElems, elems, epart);

    // Clean up
    delete[] epart;
    delete[] npart;

    return 0;
}
