#include <Eigen/Core>
#include <vector>

class Arrow {
   public:
    Arrow(const Eigen::RowVector3d& s, const Eigen::RowVector3d& e) {
        Arrow(s, e, Eigen::RowVector3d(1.0, 0, 0));
    }

    Arrow(const Eigen::RowVector3d& s, const Eigen::RowVector3d& e,
          const Eigen::RowVector3d& c)
        : start(s), end(e), color(c) {
        direction = (end - start).normalized();

        Eigen::RowVector3d per1 =
            direction.cross(Eigen::Vector3d(1, 0, 0)).normalized() * 0.5;
        if (std::isnan(per1.sum())) {
            per1 = direction.cross(Eigen::Vector3d(0, 1, 0)).normalized() * 0.5;
        }
        Eigen::RowVector3d per2 =
            direction.cross(per1.normalized()).normalized() * 0.5;

        head.resize(4);
        head[0] = end - 0.1 * (direction + per1);
        head[1] = end - 0.1 * (direction - per1);
        head[2] = end - 0.1 * (direction + per2);
        head[3] = end - 0.1 * (direction - per2);
    }

    Eigen::RowVector3d start;
    Eigen::RowVector3d end;

    Eigen::RowVector3d direction;
    std::vector<Eigen::RowVector3d> head;

    Eigen::RowVector3d color;
    size_t id;
};