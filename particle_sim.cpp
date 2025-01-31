#include <mpi.h>
#include <ascent.hpp>
#include <conduit.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <random>

// Helper to compute distance from origin
inline double distance_from_origin(double x, double y, double z)
{
    return std::sqrt(x*x + y*y + z*z);
}

int main(int argc, char* argv[])
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check args: X, R, t, u
    if(argc < 5)
    {
        if(rank == 0)
        {
            std::cerr << "Usage: mpirun -n <num_ranks> ./particle_sim "
                      << "<X> <R> <t> <u>\n"
                      << "  X = total number of particles\n"
                      << "  R = sphere radius\n"
                      << "  t = total timesteps\n"
                      << "  u = velocity magnitude\n";
        }
        MPI_Finalize();
        return 1;
    }

    // Parse arguments
    const int    X_total = std::atoi(argv[1]); // total number of particles
    const double R       = std::atof(argv[2]); // sphere radius
    const int    T       = std::atoi(argv[3]); // timesteps
    const double U       = std::atof(argv[4]); // velocity magnitude

    // Particles per rank (assuming X_total is divisible by size for simplicity)
    const int X_local = X_total / size;

    // Seed RNG differently on each rank
    std::mt19937 rng(12345 + rank);
    std::uniform_real_distribution<double> dist01(0.0, 1.0);
    std::uniform_real_distribution<double> dist_neg1_1(-1.0, 1.0);

    // Allocate storage for positions and velocities
    std::vector<double> x(X_local), y(X_local), z(X_local);
    std::vector<double> vx(X_local), vy(X_local), vz(X_local);

    // Randomly initialize particle positions inside sphere of radius R
    // We'll use a simple rejection sampling approach.
    for(int i = 0; i < X_local; i++)
    {
        while(true)
        {
            double rx = dist_neg1_1(rng);
            double ry = dist_neg1_1(rng);
            double rz = dist_neg1_1(rng);
            double rr = rx*rx + ry*ry + rz*rz;
            if(rr <= 1.0) // inside unit sphere
            {
                // scale to radius R
                x[i] = rx * R;
                y[i] = ry * R;
                z[i] = rz * R;
                break;
            }
        }
    }

    
    // Initialize random velocities with magnitude U_rand
    // random direction, and random magnitude
    for(int i = 0; i < X_local; i++)
    {
        // pick a random velocity magnitude
        const double U_rand = U * dist01(rng);
        // pick a random direction on the unit sphere
        double theta = 2.0 * M_PI * dist01(rng);
        double phi  = std::acos(1.0 - 2.0 * dist01(rng)); // 0..pi
        double sx = std::sin(phi)*std::cos(theta);
        double sy = std::sin(phi)*std::sin(theta);
        double sz = std::cos(phi);
        // scale by velocity magnitude U
        vx[i] = U_rand * sx;
        vy[i] = U_rand * sy;
        vz[i] = U_rand * sz;
    }

    // Setup Ascent
    ascent::Ascent ascent;
    conduit::Node ascent_opts;
    ascent_opts["mpi_comm"] = MPI_Comm_c2f(MPI_COMM_WORLD);
    ascent.open(ascent_opts);

    // Time integration loop
    for(int step = 0; step < T; step++)
    {
        // Update particle positions
        for(int i = 0; i < X_local; i++)
        {
            x[i] += vx[i];
            y[i] += vy[i];
            z[i] += vz[i];

            // check if out of sphere
            double r_current = distance_from_origin(x[i], y[i], z[i]);
            if(r_current > R)
            {
                // reflect the velocity by reversing normal component
                // normal direction = (x, y, z) / r_current
                double nx = x[i] / r_current;
                double ny = y[i] / r_current;
                double nz = z[i] / r_current;

                // velocity dot normal
                double vdotn = vx[i]*nx + vy[i]*ny + vz[i]*nz;
                // reflect
                vx[i] -= 2.0 * vdotn * nx;
                vy[i] -= 2.0 * vdotn * ny;
                vz[i] -= 2.0 * vdotn * nz;

                // Move back inside sphere by a small push or place on boundary
                double overshoot = r_current - R;
                // move particle back along the normal so it's exactly on R
                x[i] -= nx * overshoot;
                y[i] -= ny * overshoot;
                z[i] -= nz * overshoot;
            }
        }

        // Create a Conduit Node for our mesh data
        conduit::Node data;
        // We need to define "coordsets", "topologies", then "fields"
        // We'll store rank data under a unique domain id for blueprint
        // One approach is to set domain_id = rank.

        // Since we have point-based data:
        data["state/domain_id"] = rank;

        // Coordinates
        data["coordsets/coords/type"] = "explicit";
        data["coordsets/coords/values/x"].set_external(x);
        data["coordsets/coords/values/y"].set_external(y);
        data["coordsets/coords/values/z"].set_external(z);

        // Topology (unstructured points)
        data["topologies/particles/type"] = "points";
        data["topologies/particles/coordset"] = "coords";

        // We also want to add velocity as a vector field
        // Let's store each component
        data["fields/velocity/association"] = "vertex"; // each particle
        data["fields/velocity/topology"] = "particles";
        data["fields/velocity/volume_dependent"] = "false";
        data["fields/velocity/values/x"].set_external(vx);
        data["fields/velocity/values/y"].set_external(vy);
        data["fields/velocity/values/z"].set_external(vz);

        // For visualization, define a velocity magnitude
        {
            std::vector<double> vel_mag(X_local);
            for(int i = 0; i < X_local; i++)
            {
                vel_mag[i] = std::sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
            }
            
            data["fields/velocity_magnitude/association"] = "vertex";
            data["fields/velocity_magnitude/topology"]    = "particles";
            data["fields/velocity_magnitude/values"].set_external(vel_mag);
        }

        // Publish data to Ascent
        ascent.publish(data);

        // Execute
        conduit::Node actions;
        ascent.execute(actions);

        if(rank == 0)
        {
            std::cout << "[Step " << step << "] done\n";
        }
    }

    // Close Ascent
    ascent.close();

    MPI_Finalize();
    return 0;
}

