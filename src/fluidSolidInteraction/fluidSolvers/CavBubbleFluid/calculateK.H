const surfaceVectorField& Sf = mesh.Sf();

// Cell gradient of alpha
const tmp<volVectorField> gradAlpha = fvc::grad(alpha1);

// Interpolated face-gradient of alpha
surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);

// Face unit interface normal
surfaceVectorField nHatfv = gradAlphaf/(mag(gradAlphaf) + deltaN_);

// Face unit interface normal flux
nHatf_ = nHatfv & Sf;

// Simple expression for curvature
K_ = -fvc::div(nHatf_);
