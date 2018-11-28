#include "RigidObject.h"

RigidObject::RigidObject(const std::string& mesh_file, const ObjType t) {
    findAndLoadMesh(mesh_file);
    setType(t);
    setMass(1.0);

    // inertia of cube
    setInertia(getMass() * 2.0 / 6.0 * Eigen::Matrix3d::Identity());
    reset();
}

void RigidObject::resetMembers() {
    setLinearMomentum(Eigen::Vector3d::Zero());
    setAngularMomentum(Eigen::Vector3d::Zero());
    resetForce();
    resetTorque();
}

void RigidObject::applyForceToCOM(const Eigen::Vector3d& f) {
    if (m_type != ObjType::DYNAMIC) return;

    m_force += f;
}

void RigidObject::applyForce(const Eigen::Vector3d& f,
                             const Eigen::Vector3d& p) {
    if (m_type != ObjType::DYNAMIC) return;

    m_force += f;
    m_torque += (p - m_position).cross(f);
}

void RigidObject::applyTorque(const Eigen::Vector3d& t) {
    if (m_type != ObjType::DYNAMIC) return;

    m_torque += t;
}

void RigidObject::printDebug(const std::string& message) const {
    std::cout << std::endl;
    if (message.size() > 0) {
        std::cout << message << std::endl;
    }
    std::cout << "position:         " << m_position.transpose() << std::endl;
    std::cout << "rotation:         " << std::endl << m_rot << std::endl;
    std::cout << "linear velocity:  " << m_v.transpose() << std::endl;
    std::cout << "angular velocity: " << m_w.transpose() << std::endl;
    std::cout << "force:            " << m_force.transpose() << std::endl;
    std::cout << "torque:           " << m_torque.transpose() << std::endl;
    std::cout << "mass (inv):           " << m_mass << " (" << m_massInv << ")"
              << std::endl;
}

#pragma region GettersAndSetters
void RigidObject::setType(ObjType t) {
    m_type = t;

    if (m_type == ObjType::STATIC) {
        m_mass = std::numeric_limits<double>::infinity();
        m_massInv = 0.0;
        m_inertia.setZero();
        m_inertiaInv.setZero();
        m_force.setZero();
        m_torque.setZero();
    }
}

void RigidObject::setMass(double m) {
    if (m_type != ObjType::DYNAMIC) return;
    m_mass = m;
    m_massInv = 1.0 / m_mass;
}

void RigidObject::setInertia(const Eigen::Matrix3d& I) {
    if (m_type != ObjType::DYNAMIC) return;

    m_inertia = I;
    m_inertiaInv = m_inertia.inverse();
}

void RigidObject::setLinearMomentum(const Eigen::Vector3d& p) {
    if (m_type != ObjType::DYNAMIC) return;

    m_v = m_massInv * p;
}

void RigidObject::setAngularMomentum(const Eigen::Vector3d& l) {
    if (m_type != ObjType::DYNAMIC) return;

    m_w = getInertiaInvWorld() * l;
}

void RigidObject::setLinearVelocity(const Eigen::Vector3d& v) {
    if (m_type != ObjType::DYNAMIC) return;

    m_v = v;
}

void RigidObject::setAngularVelocity(const Eigen::Vector3d& w) {
    if (m_type != ObjType::DYNAMIC) return;

    m_w = w;
}

void RigidObject::setForce(const Eigen::Vector3d& f) {
    if (m_type != ObjType::DYNAMIC) return;

    m_force = f;
}

void RigidObject::setTorque(const Eigen::Vector3d& t) {
    if (m_type != ObjType::DYNAMIC) return;

    m_torque = t;
}

void RigidObject::resetForce() { m_force.setZero(); }

void RigidObject::resetTorque() { m_torque.setZero(); }

double RigidObject::getMass() const { return m_mass; }

double RigidObject::getMassInv() const { return m_massInv; }

Eigen::Matrix3d RigidObject::getInertia() const { return m_inertia; }

Eigen::Matrix3d RigidObject::getInertiaInv() const { return m_inertiaInv; }

Eigen::Matrix3d RigidObject::getInertiaInvWorld() const {
    return m_quat * m_inertiaInv * m_quat.inverse();
}

Eigen::Matrix3d RigidObject::getInertiaWorld() const {
    return m_quat * m_inertia * m_quat.inverse();
}

Eigen::Vector3d RigidObject::getLinearMomentum() const { return m_v * m_mass; }

Eigen::Vector3d RigidObject::getAngularMomentum() const {
    return getInertiaWorld() * m_w;
}

Eigen::Vector3d RigidObject::getLinearVelocity() const { return m_v; }

Eigen::Vector3d RigidObject::getVelocity(const Eigen::Vector3d& point) const {
    return getLinearVelocity() +
           getAngularVelocity().cross(point - getPosition());
}

Eigen::Vector3d RigidObject::getAngularVelocity() const { return m_w; }

Eigen::Vector3d RigidObject::getForce() const { return m_force; }

Eigen::Vector3d RigidObject::getTorque() const { return m_torque; }

#pragma endregion GettersAndSetters
