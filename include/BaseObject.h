#ifndef BASEOBJECT_H
#define BASEOBJECT_H

#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <Eigen/Core>

struct Mesh {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd C;

    // Per face attributes
    Eigen::MatrixXd F_normals;  // One normal per face

    // Per vertex attributes
    Eigen::MatrixXd V_normals;  // One normal per vertex

    // UV parametrization
    Eigen::MatrixXd V_uv;  // UV vertices
    Eigen::MatrixXi F_uv;  // optional faces for UVs
};

enum class ObjType { STATIC, DYNAMIC };

class BaseObject {
   public:
    bool loadMesh(const std::string& path);
    void setMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

    void findAndLoadMesh(const std::string& file);

    void reset();
    void recomputeCOM();

    void setScale(double s);
    void setID(int id);
    virtual void setType(ObjType t);
    void setPosition(const Eigen::Vector3d& p);
    void setRotation(const Eigen::Quaterniond& q);
    void setRotation(const Eigen::Matrix3d& R);
    void setColors(const Eigen::MatrixXd& C);

    double getScale() const;
    int getID() const;
    ObjType getType() const;
    Eigen::Vector3d getPosition() const;
    Eigen::Quaterniond getRotation() const;
    Eigen::Matrix3d getRotationMatrix() const;
    Eigen::Vector3d getVertexPosition(int vertexIndex) const;
    void getMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F) const;
    void getColors(Eigen::MatrixXd& C) const;

   protected:
    /*
     * Reset class variables specific to a certain object. Is called by
     * BaseObject::reset().
     */
    virtual void resetMembers() = 0;

    int m_id = -1;
    Mesh m_mesh;
    ObjType m_type;

    double m_scale = 1.0;        // Scale
    Eigen::Vector3d m_position;  // Position of the center of mass
    Eigen::Quaterniond m_quat;   // Rotation (quaternion)
    Eigen::Matrix3d m_rot;       // Rotation (matrix)
};

#endif