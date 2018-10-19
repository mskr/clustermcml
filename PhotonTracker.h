/**
* Holds photon information
* across multiple GPU runs.
*/
ALIGN_4BYTE(
struct PhotonTracker {
	float x; float y; float z; // pos [cm]
	float dx; float dy; float dz; // dir
	float weight; // 1 at start, zero when terminated
	int layerIndex; // current layer
	unsigned int rngState;
	unsigned int isDead;
});