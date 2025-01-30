#include "Structure.h"
#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <sstream>

Structure::Structure(Config *config)
{
    centerOfRotation[0] = 0.0;
    centerOfRotation[1] = 0.0;
    centerOfRotation[2] = 0.0;

    centerOfRotation_n[0] = 0.0;
    centerOfRotation_n[1] = 0.0;
    centerOfRotation_n[2] = 0.0;

    if (config->GetStructType() == "SPRING_HOR" ||
        config->GetStructType() == "SPRING_VER")
    {
        nDof = 1;
        m = config->GetSpringMass();
        Kh = config->GetSpringStiffness();
        Ch = config->GetSpringDamping();
        Ka = 0;
        Ca = 0;
        xf = 0;
        xCG = 0;
        c = 0;
        b = 0;
        S = 0;
        ICG = 0;
        If = 0;
        std::cout << "Setting mass-spring-damper system" << std::endl;
        std::cout << "Number of DOF : " << nDof << std::endl;
        std::cout << "Plunging mass : " << m << " [kg]" << std::endl;
        std::cout << "Plunging damping : " << Ch << " [Ns/m]" << std::endl;
        std::cout << "Plunging stiffness : " << Kh << " [N/m]" << std::endl;
    }
    else if (config->GetStructType() == "AIRFOIL")
    {
        nDof = 2;
        m = config->GetSpringMass();
        Kh = config->GetSpringStiffness();
        Ch = config->GetSpringDamping();
        Ka = config->GetTorsionalStiffness();
        Ca = config->GetTorsionalDamping();
        xf = config->GetFlexuralAxis();
        xCG = config->GetGravityCenter();
        c = config->GetCord();
        b = c / 2.0;
        S = m * (xCG - xf);
        // If = ICG + m*pow((xCG-xf),2);
        If = config->GetInertiaFlexural();
        // S = m*(c/2.0-xf);
        // Ia = 1.0/3.0*m*(c*c-3*c*xf+3*xf*xf);
        std::cout << "Setting pitching-plunging airfoil system" << std::endl;
        std::cout << "Number of DOF : " << nDof << std::endl;
        std::cout << "Airfoil mass : " << m << " [kg]" << std::endl;
        std::cout << "Airfoil cord : " << c << " [m]" << std::endl;
        std::cout << "Position of the flexural axis : " << xf << " [m]" << std::endl;
        std::cout << "Inertia around the flexural axis : " << If << " [kg m²]" << std::endl;
        std::cout << "Plunging damping : " << Ch << " [Ns/m]" << std::endl;
        std::cout << "Plunging stiffness : " << Kh << " [N/m]" << std::endl;
        std::cout << "Pitching damping : " << Ca << " [Ns]" << std::endl;
        std::cout << "Pitching stiffness : " << Ka << " [N]" << std::endl;
        std::cout << "Position of the center of gravity : " << xCG << " [m]" << std::endl;
        std::cout << "Static unbalance : " << S << " [kg m]" << std::endl;

        centerOfRotation[0] = xf;
        centerOfRotation_n[0] = xf;
    }
    else
    {
        nDof = 0;
        std::cerr << "Invalid structural type. Available choices are : SPRIN_HOR, SPRING_VER and AIRFOIL." << std::endl;
        throw(-1);
    }
}

void Structure::SetCenterOfRotation_X(double coord_x)
{
    centerOfRotation[0] = coord_x;
}

void Structure::SetCenterOfRotation_Y(double coord_y)
{
    centerOfRotation[1] = coord_y;
}

void Structure::SetCenterOfRotation_Z(double coord_z)
{
    centerOfRotation[2] = coord_z;
}

double Structure::GetCenterOfRotation_x() const
{
    return centerOfRotation[0];
}

double Structure::GetCenterOfRotation_y() const
{
    return centerOfRotation[1];
}

double Structure::GetCenterOfRotation_z() const
{
    return centerOfRotation[2];
}

const double *Structure::GetCenterOfRotation() const
{
    return centerOfRotation;
}

void Structure::SetCenterOfRotation_n_X(double coord_x)
{
    centerOfRotation_n[0] = coord_x;
}

void Structure::SetCenterOfRotation_n_Y(double coord_y)
{
    centerOfRotation_n[1] = coord_y;
}

void Structure::SetCenterOfRotation_n_Z(double coord_z)
{
    centerOfRotation_n[2] = coord_z;
}

double Structure::GetCenterOfRotation_n_x() const
{
    return centerOfRotation_n[0];
}

double Structure::GetCenterOfRotation_n_y() const
{
    return centerOfRotation_n[1];
}

double Structure::GetCenterOfRotation_n_z() const
{
    return centerOfRotation_n[2];
}

unsigned int Structure::GetnDof() const
{
    return nDof;
}

double Structure::Get_m() const
{
    return m;
}

double Structure::Get_Kh() const
{
    return Kh;
}

double Structure::Get_Ka() const
{
    return Ka;
}

double Structure::Get_Ch() const
{
    return Ch;
}

double Structure::Get_Ca() const
{
    return Ca;
}

double Structure::Get_S() const
{
    return S;
}

double Structure::Get_If() const
{
    return If;
}
