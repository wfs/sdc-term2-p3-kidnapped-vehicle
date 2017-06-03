/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

/** Step 1. Initialisation, with input of Initial Estimate (GPS) and Sensor Noise : see T2.L14.V3 notes for details.
 *
 * TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
 * x, y, theta and their uncertainties from GPS) and all weights to 1.
 * Add random Gaussian noise to each particle.
 *
 * NOTE: Consult particle_filter.h for more information about this method (and others in this file).
 * set number of particles
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {

    num_particles = 100;

    // define random engine and distributions used to add noise
    default_random_engine gen;
    normal_distribution<double> N_x_init(x, std[0]);
    normal_distribution<double> N_y_init(y, std[1]);
    normal_distribution<double> N_theta_init(theta, std[2]);

    for (int i = 0; i < num_particles; i++) {
        Particle p;
        p.id = i;
        p.x = N_x_init(gen);
        p.y = N_y_init(gen);
        p.theta = N_theta_init(gen);
        p.weight = 1.0; // Initialize all particles weights to 1.0
        particles.push_back(p);
    }

    is_initialized = true;
    //cout << "Init particle x: " << particles[0].x << " y: " << particles[0].y << " w: " << particles[0].weight << endl;

}

/** Step 2. Prediction, with input of Yaw Rate, Velocity and Sensor Noise. See T2.L14.V6 notes for diagram.
 *
 * From T2.L14.V6 ...
 * "Use bicycle motion model to predict where the car will be at next time step.
 * For each particle you will have to update the particle's location based on velocity and yaw rate
 * measurements to account for the uncertainty in the control input.Add gaussian noise to prediction for the project.".
 *
 * TODO: Add measurements to each particle and add random Gaussian noise.
 *
 * NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
 * http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
 * http://www.cplusplus.com/reference/random/default_random_engine/
 */
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    /*
     * Use motion model to predict where the car will be at next time step.
     * For each particle update particles location based on velocity and yaw rate
     * measurements.
     * Add gaussian noise to prediction
     */

    default_random_engine gen;

    for (int i = 0; i < num_particles; i++) {

        double x = particles[i].x;
        double y = particles[i].y;
        double theta = particles[i].theta;

        if (fabs(yaw_rate) < 1e-6) {
            x = x + velocity * delta_t * cos(theta);
            y = y + velocity * delta_t * sin(theta);

        } else {
            x = x + (velocity / yaw_rate) * (sin(theta + yaw_rate * delta_t) - sin(theta));
            y = y + (velocity / yaw_rate) * (cos(theta) - cos(theta + yaw_rate * delta_t));
            theta = theta + yaw_rate * delta_t;
        }

        // create Gaussian dist
        normal_distribution<double> N_x(x, std_pos[0] * std_pos[0]);
        normal_distribution<double> N_y(y, std_pos[1] * std_pos[1]);
        normal_distribution<double> N_yaw(theta, std_pos[2] * std_pos[2]);

        // Add random noise
        particles[i].x = N_x(gen);
        particles[i].y = N_y(gen);
        particles[i].theta = N_yaw(gen);
    }

}

/**
 * Step 3.1 Data Association. See T2.L14.V13 notes.
 *
 * TODO: Find the predicted measurement that is closest to each observed measurement and assign the
 * observed measurement to this particular landmark.
 *
 * NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
 * implement this method and use it as a helper during the updateWeights phase.
 *
 * Find the predicted measurement that is closest to each observed measurement and assign the
 * landmark Id to the observed measurement.
 */
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {

    for (int i = 0; i < observations.size(); i++) {
        double best_dist = 99999;
        double best_id = -1;

        for (int j = 0; j < predicted.size(); j++) {
            double x = observations[i].x - predicted[j].x;
            double y = observations[i].y - predicted[j].y;
            double dist = sqrt(x * x + y * y);
            if (dist < best_dist) {
                best_dist = dist;
                best_id = predicted[j].id;
            }

            // if the best distance between the predicted and observed measurements is greater that 1 m
            // do not assign to it the landmark id.
            if (best_dist <= 1.0)
                observations[i].id = best_id;
            else observations[i].id = -1;
        }
    }
}

/**
 * Step 3.0 Transformation. Update weights for each particle using a mult-variate Gaussian distribution. You can read
 * more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
 *
 * aka incorporate sensor measurements (of landmarks), such as readings from radar or lidar. See T2.L14.V11 notes.
 *
 * TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
 * more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution. See T2.L14.V13.e
 * notes for details.
 *
 * NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
 * according to the MAP'S coordinate system. You will need to transform between the two systems.
 * Keep in mind that this transformation requires both rotation AND translation (but no scaling).
 *
 * The following is a good resource for the theory:
 * https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
 * and the following is a good resource for the actual equation to implement (look at equation
 * 3.33
 * http://planning.cs.uiuc.edu/node99.html
 *
 * See Homogeneous Transformations trigonometry from my T2.L14.V13 notes links :
 * 1. Theory - http://planning.cs.uiuc.edu/node99.html
 * 2. Awesome video explanation - https://www.youtube.com/watch?v=NsiJNvsuO3s
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {

    default_random_engine gen;
    normal_distribution<double> N_x(0, std_landmark[0]);
    normal_distribution<double> N_y(0, std_landmark[1]);

    // weight normalization term
    double weight_norm = 0.0;

    double std_x = std_landmark[0];
    double std_xx = std_x * std_x;
    double std_y = std_landmark[1];
    double std_yy = std_y * std_y;
    double pnorm = 1.0 / (2.0 * M_PI * std_x * std_y);

    // iterate over each particle to find its new weight
    for (vector<Particle>::iterator it_prt = particles.begin(); it_prt != particles.end(); it_prt++) {
        double x = (*it_prt).x;
        double y = (*it_prt).y;
        double cos_theta = cos((*it_prt).theta);
        double sin_theta = sin((*it_prt).theta);

        // vector to store predicted landmark measurements
        vector<LandmarkObs> predicted;

        for (vector<Map::single_landmark_s>::iterator it_landmark = map_landmarks.landmark_list.begin();
             it_landmark != map_landmarks.landmark_list.end(); it_landmark++) {

            LandmarkObs meas;

            //Set values
            meas.id = (*it_landmark).id_i;

            // convert landmark coordinates to car coordinates, and add noise to measurement
            meas.x = (cos_theta * (*it_landmark).x_f) + (sin_theta * (*it_landmark).y_f) - (cos_theta * x) -
                     (sin_theta * y) + N_x(gen);
            meas.y = (-sin_theta * (*it_landmark).x_f) + (cos_theta * (*it_landmark).y_f) + (sin_theta * x) -
                     (cos_theta * y) + N_y(gen);

            double dist = sqrt((meas.x * meas.x) + (meas.y * meas.y));

            // select landmarks between range of measurement
            if (dist <= sensor_range) {
                predicted.push_back(meas);
            }

        }

        //associate observations to landmarks ids.
        dataAssociation(predicted, observations);

        double new_weight = 1.0;
        double prob;

        for (int i = 0; i < observations.size(); i++) {
            if (observations[i].id != -1) {

                double obsx = cos_theta * observations[i].x - sin_theta * observations[i].y + x;
                double obsy = sin_theta * observations[i].x + cos_theta * observations[i].y + y;
                double xu = obsx - map_landmarks.landmark_list[observations[i].id - 1].x_f;
                double yu = obsy - map_landmarks.landmark_list[observations[i].id - 1].y_f;

                // multi-variate gaussian distribution
                prob = pnorm * exp(-((xu * xu) / (2 * std_xx) + (yu * yu) / (2 * std_yy)));
                new_weight *= prob;

            } else {
                new_weight = 0.0;
                break;
            }
        }

        (*it_prt).weight = new_weight;
        weight_norm += new_weight;

    }

    // normalise new particles weights
    for (vector<Particle>::iterator it_prt = particles.begin(); it_prt != particles.end(); it_prt++) {
        (*it_prt).weight = (*it_prt).weight / weight_norm;

    }

}

/**
 * Step 4. Resample the particles with the probability proportional to these weights. See T2.L14.V11 notes.
 *
 * TODO: Resample particles with replacement with probability proportional to their weight.
 *
 * NOTE: You may find std::discrete_distribution helpful here.
 * http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
 *
 */
void ParticleFilter::resample() {

    // Set of current particles
    std::vector<Particle> particles_tmp;

    std::default_random_engine int_eng{};
    std::uniform_int_distribution<> int_distribution{0, num_particles}; // type of engine
    std::default_random_engine real_eng{};
    std::uniform_real_distribution<> real_distribution{0, 1}; // type of engine

    int index = int_distribution(int_eng);

    double beta = 0.0;

    double nw = 0;

    for (int i = 0; i < particles.size(); i++) {
        if (nw < particles[i].weight)
            nw = particles[i].weight;
    }

    // See Sebastian Thrun's T2.L13.V20 Resampling Wheel algorithm explanation and
    // T2.L13.V20.nn for pseudo code correction details.
    for (int i = 0; i < num_particles; i++) {
        beta += real_distribution(real_eng) * 2.0 * nw;
        while (beta > particles[index].weight) {
            beta -= particles[index].weight;
            index = (index + 1) % num_particles;
        }
        particles_tmp.push_back(particles[index]);

    }
    // set resampled particles
    particles = particles_tmp;

}

Particle
ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x,
                                std::vector<double> sense_y) {
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best) {
    vector<int> v = best.associations;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best) {
    vector<double> v = best.sense_x;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best) {
    vector<double> v = best.sense_y;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}
