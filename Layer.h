#include "clmem.h"

ALIGN_4BYTE(
struct Layer {
	float absorbCoeff;
	float scatterCoeff;
	float g; // anisotropy
	float n; // refractive index
});