#include "BaseObject.h"
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>

bool BaseObject::loadMesh(const std::string& path) {
    bool succ = false;

    std::ifstream infile(path);
    if (!infile.good()) {
        return false;
    }

    const std::string OFF(".off");
    if (path.compare(path.size() - OFF.size(), OFF.size(), OFF) == 0) {
        succ = igl::readOFF(path, m_mesh.V, m_mesh.F, m_mesh.V_normals);
        if (succ) {
            std::cout << "Reading OFF-file from " << path << " ..."
                      << std::endl;
        }
    }

    const std::string OBJ(".obj");
    if (path.compare(path.size() - OBJ.size(), OBJ.size(), OBJ) == 0) {
        succ = igl::readOBJ(path, m_mesh.V, m_mesh.F);
        if (succ) {
            std::cout << "Reading OBJ-file from " << path << " ..."
                      << std::endl;
            igl::per_vertex_normals(m_mesh.V, m_mesh.F, m_mesh.V_normals);
        }
    }
    m_mesh.C = Eigen::MatrixXd(1, 3);
    m_mesh.C << 255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0;
    return succ;
}

void BaseObject::setMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    m_mesh.V = V;
    m_mesh.F = F;
}

void BaseObject::findAndLoadMesh(const std::string& file) {
    if (loadMesh(file)) return;
    if (loadMesh("data/" + file)) return;
    if (loadMesh("../data/" + file)) return;
    if (loadMesh("../../data/" + file)) return;
    if (loadMesh("../../../data/" + file)) return;
    std::cerr << "Failed to find " << file << std::endl;
}

void BaseObject::reset() {
    setPosition(Eigen::Vector3d::Zero());
    setRotation(Eigen::Matrix3d::Identity());
    resetMembers();
}

void BaseObject::recomputeCOM() {
    Eigen::Vector3d COM = m_mesh.V.colwise().mean();
    m_mesh.V = m_mesh.V.rowwise() - COM.transpose();
}

void BaseObject::setScale(double s) { m_scale = s; }

void BaseObject::setID(int id) { m_id = id; }

void BaseObject::setType(ObjType t) { m_type = t; }

void BaseObject::setPosition(const Eigen::Vector3d& p) { m_position = p; }

void BaseObject::setRotation(const Eigen::Quaterniond& q) {
    m_quat = q;
    m_rot = q.toRotationMatrix();
}

void BaseObject::setRotation(const Eigen::Matrix3d& R) {
    m_rot = R;
    m_quat = R;
}

void BaseObject::setColors(const Eigen::MatrixXd& C) { m_mesh.C = C; }

double BaseObject::getScale() const { return m_scale; }

int BaseObject::getID() const { return m_id; }

ObjType BaseObject::getType() const { return m_type; }

Eigen::Vector3d BaseObject::getPosition() const { return m_position; }

Eigen::Quaterniond BaseObject::getRotation() const { return m_quat; }

Eigen::Matrix3d BaseObject::getRotationMatrix() const { return m_rot; }

Eigen::Vector3d BaseObject::getVertexPosition(int vertexIndex) const {
    return m_mesh.V.row(vertexIndex) * m_scale *
               getRotationMatrix().transpose() +
           getPosition().transpose();
}

void BaseObject::getMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F) const {
    // get mesh after rotation and translation
    V = (m_mesh.V * m_scale * getRotationMatrix().transpose()).rowwise() +
        getPosition().transpose();
    F = m_mesh.F;
}

void BaseObject::getColors(Eigen::MatrixXd& C) const { C = m_mesh.C; }