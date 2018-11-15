#include "CannonBallSim.h"

/////////////////////////////////////
//////// EX 1 - PROBLEM 1 ///////////
/////////////////////////////////////
bool CannonBallSim::advance() {
    // perform time integration with different integrators

	// use p_ball, m_dt, m_gravity
    switch (m_method) {
        case 0:
            {
            Eigen::Vector3d v = p_ball->getLinearVelocity();

            Eigen::Vector3d pos = p_ball->getPosition();
            pos[0] += v[0] * m_dt;
            pos[1] += v[1] * m_dt + .5*m_gravity[1]*m_dt*m_dt;
            pos[2] += v[2] * m_dt;
            p_ball->setPosition(pos);

            v[1] += m_gravity[1] * m_dt;
            p_ball->setLinearVelocity(v);
            }
            break;

        case 1:
            // explicit euler
            {
            Eigen::Vector3d v = p_ball->getLinearVelocity();

            Eigen::Vector3d pos = p_ball->getPosition();
            pos[0] += v[0] * m_dt;
            pos[1] += v[1] * m_dt;
            pos[2] += v[2] * m_dt;
            p_ball->setPosition(pos);

            v[1] += m_gravity[1] * m_dt;
            p_ball->setLinearVelocity(v);
            }
            break;

        case 2:
            // symplectic euler
            {
            Eigen::Vector3d v = p_ball->getLinearVelocity();
            v[1] += m_gravity[1] * m_dt;
            p_ball->setLinearVelocity(v);

            Eigen::Vector3d pos = p_ball->getPosition();
            pos[0] += v[0] * m_dt;
            pos[1] += v[1] * m_dt;
            pos[2] += v[2] * m_dt;
            p_ball->setPosition(pos);
            }
            break;

        default:
            std::cerr << m_method << " is not a valid integrator method."
                        << std::endl;
    }

    // advance time
    m_time += m_dt;
    m_step++;

    // log
    if ((m_step % m_log_frequency) == 0) {
        m_trajectories.back().push_back(p_ball->getPosition());
    }

    return false;
}
