/** This program demonstrates how
 *  to partition an unstructured 
 *  mesh using METIS.
 *
 *  It partitions a mesh using its dual graph.
 *
 *  @author: Andrew C. Kirby
 *  @data: July 1, 2024
 *  @info: IHPCSS 2024, Kobe, Japan
 */

/* system header files */
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

/* library header files */
#include <metis.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

// Struct for storing node coordinates
struct Node {
    double x, y;
};

// Struct for storing an element (e.g., triangle) defined by its nodes
struct Element {
    int nodes[3]; // For triangles
};

// Struct for storing an unstructured mesh
struct Mesh {
    idx_t numNodes;
    idx_t numElems;
    idx_t numEdges;
    std::vector<Node> nodes;
    std::vector<Element> elems;
    std::vector<idx_t> npart;
    std::vector<idx_t> epart;
};

void readMesh(const std::string filename,Mesh &mesh){
    std::ifstream file(filename);

    // Skip 1st line comments
    std::string comments;
    std::getline(file, comments);

    // Read node, element, and edge counts
    file >> mesh.numNodes >> mesh.numElems >> mesh.numEdges;
    std::cout << ANSI_COLOR_GREEN 
                 "MESH STATS:\n" ANSI_COLOR_RESET
              << "  # Nodes : " << mesh.numNodes << "\n"
              << "  # Elems : " << mesh.numElems << "\n"
              << "  # Edges : " << mesh.numEdges << std::endl; 

    // Allocate node and element lists
    mesh.elems.resize(mesh.numElems);
    mesh.nodes.resize(mesh.numNodes);

    // Read node coordinates
    for (int i = 0; i < mesh.numNodes; ++i) {
        file >> mesh.nodes[i].x >> mesh.nodes[i].y;
    }

    int n0,n1,n2;
    // Read element node IDs
    for (int i = 0; i < mesh.numElems; ++i) {
        file >> n0 >> n1 >> n2;

        // make nodes base 0
        mesh.elems[i].nodes[0] = n0 - 1;
        mesh.elems[i].nodes[1] = n1 - 1;
        mesh.elems[i].nodes[2] = n2 - 1;
   }
}

// Function to write mesh in VTK format
void writeMeshVTK(const std::string filename,
                  int numNodes, Node* nodes,
                  int numElems, Element* elems,
                  idx_t* part) {
    std::ofstream file(filename);

    // Write VTK header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "Partitioned Mesh Example" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write nodes
    file << "POINTS " << numNodes << " double" << std::endl;
    for (int i = 0; i < numNodes; ++i) {
        file << nodes[i].x << " " << nodes[i].y << " 0.0" << std::endl; // Assuming 2D mesh (z=0.0)
    }

    // Write elements
    int totalNodes = 0;
    for (int i = 0; i < numElems; ++i) {
        totalNodes += 1 + elems[i].nodes[2]; // Number of nodes in the current element
    }
    file << "CELLS " << numElems << " " << totalNodes << std::endl;
    for (int i = 0; i < numElems; ++i) {
        file << "3 "; // Three nodes per triangle (assuming triangles)
        file << elems[i].nodes[0] << " " << elems[i].nodes[1] << " " << elems[i].nodes[2] << std::endl;
    }

    // Write element types
    file << "CELL_TYPES " << numElems << std::endl;
    for (int i = 0; i < numElems; ++i) {
        file << "5\n"; // VTK_TRIANGLE type
    }
    file << std::endl;

    // Write partition data as a scalar field
    file << "CELL_DATA " << numElems << std::endl;
    file << "SCALARS Partition_ID int 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < numElems; ++i) {
        file << part[i] << std::endl;
    }

    file.close();
    std::cout << "VTK File Written: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    Mesh mesh;
    
    int meshChoice = 1;
    int numPartitions = 2;
    std::string meshTag  = ".mesh";
    std::string vtkTag   = ".vtk";

    // read number of partitions requested
    if(argc <= 1) {
        std::cout << "USAGE" << std::endl;
        std::cout << ANSI_COLOR_GREEN  "  ./MetisDemo "
                     ANSI_COLOR_CYAN   "<MeshID{int}> "
                     ANSI_COLOR_YELLOW "<NParts{int}>\n"
                     ANSI_COLOR_RESET  "     OPTIONS:  " << ANSI_COLOR_CYAN << "MeshID: " ANSI_COLOR_RESET "[1] Box, [2] Airfoil\n"
                     ANSI_COLOR_YELLOW "               NParts:" ANSI_COLOR_RESET " >= 2\n\n";
    }
    if(argc > 1) meshChoice = std::atoi(argv[1]); 
    if(argc > 2) numPartitions = std::atoi(argv[2]);

    std::string baseName = (meshChoice==1) ? "box":"naca0012";
    std::string meshFile = baseName + meshTag;

    std::cout << "Mesh: " << meshFile << std::endl;
    std::cout << "Partitions: " << numPartitions << std::endl;
    if (numPartitions <= 1) {
        std::cout << "ERROR: Please enter >=2 partitions." << std::endl;
        exit(0);
    }

    // read mesh file
    std::string vtkFile  = baseName + ".parts" + std::to_string(numPartitions) + vtkTag;
    readMesh(meshFile,mesh);

    const int nPi = 3; 
    idx_t *eptr = new idx_t[mesh.numElems + 1];
    eptr[0] = 0;
    for(int i = 0 ; i < mesh.numElems; i++) eptr[i+1] = eptr[i] + nPi;
    
    idx_t *eind = new idx_t[eptr[mesh.numElems]];
    int cnt = 0;
    for (int i = 0 ; i < mesh.numElems; i++){
        for (int j = 0; j < nPi; j++){
            eind[cnt] = mesh.elems[i].nodes[j];
            cnt++;
        }
    }

    // Convert the mesh to METIS-compatible format
    idx_t ncommon = 2;  // Number of common nodes that two elements must have in order to an edge between them in the dual graph.
    idx_t numParts = (idx_t) numPartitions; // Number of partitions
    
    mesh.epart.resize(mesh.numElems);
    mesh.npart.resize(mesh.numNodes);

    // METIS options
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // Partition the mesh using METIS
    idx_t objval; // METIS output (total edge-cut or communication volume)
    METIS_PartMeshDual(&mesh.numElems,     // number of mesh elements
                       &mesh.numNodes,     // number of mesh nodes
                       eptr,               // array storing mesh node counts per element (scan)
                       eind,               // array storing mesh info
                       NULL,               // vwgt: An array of size ne specifying the weights of the elements.
                       NULL,               // vsize: 
                       &ncommon,           // number of common nodes two elements need to share a face
                       &numParts,          // number of partitions to generate
                       NULL,               // tpwghts: desired weight for each partition
                       options,            // METIS Options
                       &objval,            // METIS Output
                       mesh.epart.data(),  // partition number of each element
                       mesh.npart.data()); // partition number of each node
    
    // Output partitioned mesh in VTK format
    writeMeshVTK(vtkFile,
                 mesh.numNodes, mesh.nodes.data(), 
                 mesh.numElems, mesh.elems.data(),
                 mesh.epart.data());

    return 0;
}
