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
	unsigned int rngState; // keep rng state to avoid reusing seeds
	unsigned int isDead; // mark as dead when pipeline contains enough photons
});