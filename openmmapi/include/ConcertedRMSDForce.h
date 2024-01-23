#ifndef OPENMM_CONCERTEDRMSDFORCE_H_
#define OPENMM_CONCERTEDRMSDFORCE_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2024 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/Force.h"
#include "openmm/Vec3.h"
#include <vector>
#include "internal/windowsExportOpenMMLab.h"

using namespace OpenMM;

namespace OpenMMLab {

/**
 * This is a force whose energy equals the root mean squared deviation (RMSD)
 * between the current coordinates and a reference structure.  It is intended for
 * use with CustomCVForce.  You will not normally want a force that exactly equals
 * the RMSD, but there are many situations where it is useful to have a restraining
 * or biasing force that depends on the RMSD in some way.
 *
 * The force is computed by first aligning the particle positions to the reference
 * structure, then computing the RMSD between the aligned positions and the reference.
 * The computation can optionally be done based on only a subset of the particles
 * in the system.
 */

class OPENMM_EXPORT_OPENMM_LAB ConcertedRMSDForce : public Force {
public:
    /**
     * Create an ConcertedRMSDForce.
     *
     * @param referencePositions  the reference positions to compute the deviation
     *                            from.  The length of this vector must equal the
     *                            number of particles in the system, even if not
     *                            all particles are used in computing the RMSD.
     * @param particles           the indices of the particles to use when computing
     *                            the RMSD.  If this is empty (the default), all
     *                            particles in the system will be used.
     */
    explicit ConcertedRMSDForce(const std::vector<Vec3>& referencePositions,
                       const std::vector<int>& particles=std::vector<int>());
    /**
     * Get the reference positions to compute the deviation from.
     */
    const std::vector<Vec3>& getReferencePositions() const {
        return referencePositions;
    }
    /**
     * Set the reference positions to compute the deviation from.
     */
    void setReferencePositions(const std::vector<Vec3>& positions);
    /**
     * Get the indices of the particles to use when computing the RMSD.  If this
     * is empty, all particles in the system will be used.
     */
    const std::vector<int>& getParticles() const {
        return particles;
    }
    /**
     * Set the indices of the particles to use when computing the RMSD.  If this
     * is empty, all particles in the system will be used.
     */
    void setParticles(const std::vector<int>& particles);
    /**
     * Update the reference positions and particle indices in a Context to match those stored
     * in this Force object.  This method provides an efficient method to update certain parameters
     * in an existing Context without needing to reinitialize it.  Simply call setReferencePositions()
     * and setParticles() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
protected:
    ForceImpl* createImpl() const;
private:
    std::vector<Vec3> referencePositions;
    std::vector<int> particles;
};

} // namespace OpenMMLab

#endif /*OPENMM_CONCERTEDRMSDFORCE_H_*/
