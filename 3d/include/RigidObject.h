#ifndef RIGIDOBJECT_H
#define RIGIDOBJECT_H

#include "BaseObject.h"

/*
 * Base class representing a simple rigid object.
 */
class RigidObject : public BaseObject {
   public:
    RigidObject() {}
    RigidObject(const std::string& mesh_path,
                const ObjType t = ObjType::DYNAMIC);

    void applyForceToCOM(const Eigen::Vector3d& f);
    /*
     * Apply force f at point p.
     */
    void applyForce(const Eigen::Vector3d& f, const Eigen::Vector3d& p);
    void applyTorque(const Eigen::Vector3d& t);
    void printDebug(const std::string& message = "") const;

#pragma region GettersAndSetters
    virtual void setType(ObjType t);
    void setMass(double m);
    void setInertia(const Eigen::Matrix3d& I);
    void setLinearMomentum(const Eigen::Vector3d& p);
    void setAngularMomentum(const Eigen::Vector3d& l);
    void setLinearVelocity(const Eigen::Vector3d& v);
    void setAngularVelocity(const Eigen::Vector3d& w);
    void setForce(const Eigen::Vector3d& f);
    void setTorque(const Eigen::Vector3d& t);
    void resetForce();
    void resetTorque();

    double getMass() const;
    double getMassInv() const;
    Eigen::Matrix3d getInertia() const;
    Eigen::Matrix3d getInertiaInv() const;
    Eigen::Matrix3d getInertiaInvWorld() const;
    Eigen::Matrix3d getInertiaWorld() const;
    Eigen::Vector3d getLinearMomentum() const;
    Eigen::Vector3d getAngularMomentum() const;
    Eigen::Vector3d getLinearVelocity() const;
    Eigen::Vector3d getVelocity(const Eigen::Vector3d& point) const;
    Eigen::Vector3d getAngularVelocity() const;
    Eigen::Vector3d getForce() const;
    Eigen::Vector3d getTorque() const;
#pragma endregion GettersAndSetters

   protected:
    virtual void resetMembers() override;

   private:
    double m_mass;                 // Body mass
    double m_massInv;              // Inverted mass
    Eigen::Matrix3d m_inertia;     // Intertial Tensor (initially set to cube)
    Eigen::Matrix3d m_inertiaInv;  // Inverse

    Eigen::Vector3d m_v;  // Linear velocity
    Eigen::Vector3d m_w;  // Angular velocity

    Eigen::Vector3d m_force;   // Force on body
    Eigen::Vector3d m_torque;  // Torque on body
};

#endif