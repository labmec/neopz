/**
 * @file
 * @brief Contains the implementation of the ConvTest methods. 
 */
#include "convtest.h"

ConvTest::ConvTest()
{
}


ConvTest::~ConvTest()
{
}

#include "pzgeoel.h"

/**Avaliate the Jacobian by Obtained Convergence Order*/
void ConvTest::JacobianConv(TPZGeoEl &Object, TPZVec< REAL > StartPoint)
{
	TPZVec< REAL > Out(3,0.);
	if(Object.Dimension() == 1)
	{
		TPZFMatrix<REAL> jacobian(1,1);
		TPZFMatrix<REAL> Axes(1,3);
		REAL detJacobian;
		TPZFMatrix<REAL> InvJac(1,1);
		TPZVec< REAL > QsiEtaIni(1);
		QsiEtaIni[0] = StartPoint[0];
		const double deltaQsi = 0.1;
		double alpha;
		
		std::cout << "\ninitial Qsi = " << QsiEtaIni[0] << "\n";
		std::cout << "deltaQsi = const = " << deltaQsi << "\n\n";
		
		TPZVec< REAL > OutAprox(3);
		TPZFMatrix<REAL> error(11,1,0.);
		double dX, dY, dZ;
		
		Object.Jacobian(QsiEtaIni,jacobian,Axes,detJacobian,InvJac);
		
		for(int i = 0; i < 11; i++)
		{
			alpha = i/10.;
			Object.X(QsiEtaIni,Out);//for aproximate compute
			//                std::cout << "\nalpha = " << alpha << std::endl;
			//                std::cout << "f(x)                     : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			
			dX = alpha*( jacobian.GetVal(0,0)*deltaQsi )*Axes(0,0);
			OutAprox[0] = Out[0] + dX;
			
			dY = alpha*( jacobian.GetVal(0,0)*deltaQsi )*Axes(0,1);
			OutAprox[1] = Out[1] + dY;
			
			dZ = alpha*( jacobian.GetVal(0,0)*deltaQsi )*Axes(0,2);
			OutAprox[2] = Out[2] + dZ;
			
			StartPoint[0] = QsiEtaIni[0] + alpha*deltaQsi;
			
			Object.X(StartPoint,Out);//for real compute
			
			//                std::cout << "alpha.(axes.J).dx        : x = " << dX << " | y = " << dY << " | z = " << dZ << "\n";
			//                std::cout << "f(x) + alpha.(axes.J).dx : x = " << OutAprox[0] << " | y = " << OutAprox[1] << " | z = " << OutAprox[2] << "\n";
			//                std::cout << "f(x + alpha.dx)          : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			//                std::cout << "--------------------------------------------------------------------\n";
			
			Out[0] -= OutAprox[0];
			Out[1] -= OutAprox[1];
			Out[2] -= OutAprox[2];
			
			double XDiffNorm = sqrt(Out[0]*Out[0] + Out[1]*Out[1] + Out[2]*Out[2]);
			error(int(i),0) = XDiffNorm;
		}
		
		//           std::cout << "ERROR Vector:\n"; error.Print();
		
		std::cout << "Convergence Order:\n";
		for(int j = 2; j < 11; j++)
		{
			std::cout << ( log(error(1,0)) - log(error(j,0)) )/( log(0.1)-log(j/10.) ) << "\n";
		}
		std::cout << "\nIf another kind of results are needed, edit the ConvTest Class on source code!\n";
		return;
	}
	
	if(Object.Dimension() == 2)
	{
		TPZFMatrix<REAL> jacobian(2,2);
		TPZFMatrix<REAL> Axes(2,3);
		REAL detJacobian;
		TPZFMatrix<REAL> InvJac(2,2);
		TPZVec< REAL > QsiEtaIni(2);
		QsiEtaIni[0] = StartPoint[0];
		QsiEtaIni[1] = StartPoint[1];
		const double deltaQsi = 0.1;
		const double deltaEta = 0.1;
		double alpha;
		
		std::cout << "\ninitial Qsi = " << QsiEtaIni[0] << " | initial Eta = " << QsiEtaIni[1] << "\n";
		std::cout << "deltaQsi = const = " << deltaQsi << " | deltaEta = const = " << deltaEta << "\n\n";
		
		TPZVec< REAL > OutAprox(3);
		TPZFMatrix<REAL> error(11,1,0.);
		double dX, dY, dZ;
		
		Object.Jacobian(QsiEtaIni,jacobian,Axes,detJacobian,InvJac);
		
		for(int i = 0; i <= 10; i++)
		{
			alpha = i/10.;
			Object.X(QsiEtaIni,Out);//for aproximate compute
			//                std::cout << "\nalpha = " << alpha << std::endl;
			//                std::cout << "f(x)                     : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			
			dX = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta )*Axes(0,0) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta )*Axes(1,0);
			OutAprox[0] = Out[0] + dX;
			
			dY = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta )*Axes(0,1) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta )*Axes(1,1);
			OutAprox[1] = Out[1] + dY;
			
			dZ = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta )*Axes(0,2) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta )*Axes(1,2);
			OutAprox[2] = Out[2] + dZ;
			
			StartPoint[0] = QsiEtaIni[0] + alpha*deltaQsi;
			StartPoint[1] = QsiEtaIni[1] + alpha*deltaEta;
			Object.X(StartPoint,Out);//for real compute
			
			//                std::cout << "alpha.(axes.J).dx        : x = " << dX << " | y = " << dY << " | z = " << dZ << "\n";
			//                std::cout << "f(x) + alpha.(axes.J).dx : x = " << OutAprox[0] << " | y = " << OutAprox[1] << " | z = " << OutAprox[2] << "\n";
			//                std::cout << "f(x + alpha.dx)          : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			//                std::cout << "--------------------------------------------------------------------\n";
			
			Out[0] -= OutAprox[0];
			Out[1] -= OutAprox[1];
			Out[2] -= OutAprox[2];
			
			double XDiffNorm = sqrt(Out[0]*Out[0] + Out[1]*Out[1] + Out[2]*Out[2]);
			error(int(i),0) = XDiffNorm;
		}
		
		//           std::cout << "ERROR Vector:\n"; error.Print();
		
		std::cout << "Convergence Order:\n";
		for(int j = 2; j < 11; j++)
		{
			std::cout << ( log(error(1,0)) - log(error(j,0)) )/( log(0.1)-log(j/10.) ) << "\n";
		}
		std::cout << "\nIf another kind of results are needed, edit the ConvTest Class on source code!\n";
		return;
	}
	
	if(Object.Dimension() == 3)
	{
		TPZFMatrix<REAL> jacobian(3,3);
		TPZFMatrix<REAL> Axes(3,3);
		REAL detJacobian;
		TPZFMatrix<REAL> InvJac(3,3);
		TPZVec< REAL > QsiEtaIni(3);
		QsiEtaIni[0] = StartPoint[0];
		QsiEtaIni[1] = StartPoint[1];
		QsiEtaIni[2] = StartPoint[2];
		const double deltaQsi  = 0.1;
		const double deltaEta  = 0.1;
		const double deltaZeta = 0.1;
		double alpha;
		
		std::cout << "\ninitial Qsi = " << QsiEtaIni[0] << " | initial Eta = " << QsiEtaIni[1] << " | initial Zeta = " << QsiEtaIni[2] <<"\n";
		std::cout << "deltaQsi = const = " << deltaQsi << " | deltaEta = const = " << deltaEta << " | deltaZeta = const = " << deltaZeta <<"\n\n";
		
		TPZVec< REAL > OutAprox(3);
		TPZFMatrix<REAL> error(11,1,0.);
		double dX, dY, dZ;
		
		Object.Jacobian(QsiEtaIni,jacobian,Axes,detJacobian,InvJac);
		
		for(int i = 0; i <= 10; i++)
		{
			alpha = i/10.;
			Object.X(QsiEtaIni,Out);//for aproximate compute
			//               std::cout << "\nalpha = " << alpha << std::endl;
			//               std::cout << "f(x)                     : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			
			dX = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta + jacobian.GetVal(0,2)*deltaZeta)*Axes(0,0) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta + jacobian.GetVal(1,2)*deltaZeta)*Axes(1,0) + alpha*( jacobian.GetVal(2,0)*deltaQsi + jacobian.GetVal(2,1)*deltaEta + jacobian.GetVal(2,2)*deltaZeta)*Axes(2,0);
			OutAprox[0] = Out[0] + dX;
			
			dY = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta + jacobian.GetVal(0,2)*deltaZeta)*Axes(0,1) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta + jacobian.GetVal(1,2)*deltaZeta)*Axes(1,1) + alpha*( jacobian.GetVal(2,0)*deltaQsi + jacobian.GetVal(2,1)*deltaEta + jacobian.GetVal(2,2)*deltaZeta)*Axes(2,1);
			OutAprox[1] = Out[1] + dY;
			
			dZ = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta + jacobian.GetVal(0,2)*deltaZeta)*Axes(0,2) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta + jacobian.GetVal(1,2)*deltaZeta)*Axes(1,2) + alpha*( jacobian.GetVal(2,0)*deltaQsi + jacobian.GetVal(2,1)*deltaEta + jacobian.GetVal(2,2)*deltaZeta)*Axes(2,2);
			OutAprox[2] = Out[2] + dZ;
			
			StartPoint[0] = QsiEtaIni[0] + alpha*deltaQsi;
			StartPoint[1] = QsiEtaIni[1] + alpha*deltaEta;
			StartPoint[2] = QsiEtaIni[2] + alpha*deltaZeta;
			Object.X(StartPoint,Out);//for real compute
			
			//               std::cout << "alpha.(axes.J).dx        : x = " << dX << " | y = " << dY << " | z = " << dZ << "\n";
			//               std::cout << "f(x) + alpha.(axes.J).dx : x = " << OutAprox[0] << " | y = " << OutAprox[1] << " | z = " << OutAprox[2] << "\n";
			//               std::cout << "f(x + alpha.dx)          : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			//               std::cout << "--------------------------------------------------------------------\n";
			
			Out[0] -= OutAprox[0];
			Out[1] -= OutAprox[1];
			Out[2] -= OutAprox[2];
			
			double XDiffNorm = sqrt(Out[0]*Out[0] + Out[1]*Out[1] + Out[2]*Out[2]);
			error(int(i),0) = XDiffNorm;
		}
		
		//           std::cout << "ERROR Vector:\n"; error.Print();
		
		std::cout << "Convergence Order:\n";
		for(int j = 2; j < 11; j++)
		{
			std::cout << ( log(error(1,0)) - log(error(j,0)) )/( log(0.1)-log(j/10.) ) << "\n";
		}
		std::cout << "\nIf another kind of results are needed, edit the ConvTest Class on source code!\n";
		return;
	}
	else
	{
		std::cout << "Element don't fit in an option of Convergence Analys!\nSee ConvTest Class...\n";
		return;
	}
}

/**Avaliate the Jacobian by Obtained Convergence Order*/
void ConvTest::JacobianConv(TPZGeoElSide &Object, TPZVec< REAL > StartPoint)
{
	TPZVec< REAL > Out(3,0.);
	if(Object.Dimension() == 1)
	{
		TPZFMatrix<REAL> jacobian(1,1);
		TPZFMatrix<REAL> Axes(1,3);
		REAL detJacobian;
		TPZFMatrix<REAL> InvJac(1,1);
		TPZVec< REAL > QsiEtaIni(1);
		QsiEtaIni[0] = StartPoint[0];
		const double deltaQsi = 0.1;
		double alpha;
		
		std::cout << "\ninitial Qsi = " << QsiEtaIni[0] << "\n";
		std::cout << "deltaQsi = const = " << deltaQsi << "\n\n";
		
		TPZVec< REAL > OutAprox(3);
		TPZFMatrix<REAL> error(11,1,0.);
		double dX, dY, dZ;
		
		Object.Jacobian(QsiEtaIni,jacobian,Axes,detJacobian,InvJac);
		
		for(int i = 0; i < 11; i++)
		{
			alpha = i/10.;
			Object.X(QsiEtaIni,Out);//for aproximate compute
			//                std::cout << "\nalpha = " << alpha << std::endl;
			//                std::cout << "f(x)                     : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			
			dX = alpha*( jacobian.GetVal(0,0)*deltaQsi )*Axes(0,0);
			OutAprox[0] = Out[0] + dX;
			
			dY = alpha*( jacobian.GetVal(0,0)*deltaQsi )*Axes(0,1);
			OutAprox[1] = Out[1] + dY;
			
			dZ = alpha*( jacobian.GetVal(0,0)*deltaQsi )*Axes(0,2);
			OutAprox[2] = Out[2] + dZ;
			
			StartPoint[0] = QsiEtaIni[0] + alpha*deltaQsi;
			
			Object.X(StartPoint,Out);//for real compute
			
			//                std::cout << "alpha.(axes.J).dx        : x = " << dX << " | y = " << dY << " | z = " << dZ << "\n";
			//                std::cout << "f(x) + alpha.(axes.J).dx : x = " << OutAprox[0] << " | y = " << OutAprox[1] << " | z = " << OutAprox[2] << "\n";
			//                std::cout << "f(x + alpha.dx)          : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			//                std::cout << "--------------------------------------------------------------------\n";
			
			Out[0] -= OutAprox[0];
			Out[1] -= OutAprox[1];
			Out[2] -= OutAprox[2];
			
			double XDiffNorm = sqrt(Out[0]*Out[0] + Out[1]*Out[1] + Out[2]*Out[2]);
			error(int(i),0) = XDiffNorm;
		}
		
		//           std::cout << "ERROR Vector:\n"; error.Print();
		
		std::cout << "Convergence Order:\n";
		for(int j = 2; j < 11; j++)
		{
			std::cout << ( log(error(1,0)) - log(error(j,0)) )/( log(0.1)-log(j/10.) ) << "\n";
		}
		std::cout << "\nIf another kind of results are needed, edit the ConvTest Class on source code!\n";
		return;
	}
	
	if(Object.Dimension() == 2)
	{
		TPZFMatrix<REAL> jacobian(2,2);
		TPZFMatrix<REAL> Axes(2,3);
		REAL detJacobian;
		TPZFMatrix<REAL> InvJac(2,2);
		TPZVec< REAL > QsiEtaIni(2);
		QsiEtaIni[0] = StartPoint[0];
		QsiEtaIni[1] = StartPoint[1];
		const double deltaQsi = 0.1;
		const double deltaEta = 0.1;
		double alpha;
		
		std::cout << "\ninitial Qsi = " << QsiEtaIni[0] << " | initial Eta = " << QsiEtaIni[1] << "\n";
		std::cout << "deltaQsi = const = " << deltaQsi << " | deltaEta = const = " << deltaEta << "\n\n";
		
		TPZVec< REAL > OutAprox(3);
		TPZFMatrix<REAL> error(11,1,0.);
		double dX, dY, dZ;
		
		Object.Jacobian(QsiEtaIni,jacobian,Axes,detJacobian,InvJac);
		
		for(int i = 0; i <= 10; i++)
		{
			alpha = i/10.;
			Object.X(QsiEtaIni,Out);//for aproximate compute
			//                std::cout << "\nalpha = " << alpha << std::endl;
			//                std::cout << "f(x)                     : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			
			dX = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta )*Axes(0,0) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta )*Axes(1,0);
			OutAprox[0] = Out[0] + dX;
			
			dY = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta )*Axes(0,1) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta )*Axes(1,1);
			OutAprox[1] = Out[1] + dY;
			
			dZ = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta )*Axes(0,2) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta )*Axes(1,2);
			OutAprox[2] = Out[2] + dZ;
			
			StartPoint[0] = QsiEtaIni[0] + alpha*deltaQsi;
			StartPoint[1] = QsiEtaIni[1] + alpha*deltaEta;
			Object.X(StartPoint,Out);//for real compute
			
			//                std::cout << "alpha.(axes.J).dx        : x = " << dX << " | y = " << dY << " | z = " << dZ << "\n";
			//                std::cout << "f(x) + alpha.(axes.J).dx : x = " << OutAprox[0] << " | y = " << OutAprox[1] << " | z = " << OutAprox[2] << "\n";
			//                std::cout << "f(x + alpha.dx)          : x = " << Out[0] << " | y = " << Out[1] << " | z = " << Out[2] << "\n";
			//                std::cout << "--------------------------------------------------------------------\n";
			
			Out[0] -= OutAprox[0];
			Out[1] -= OutAprox[1];
			Out[2] -= OutAprox[2];
			
			double XDiffNorm = sqrt(Out[0]*Out[0] + Out[1]*Out[1] + Out[2]*Out[2]);
			error(int(i),0) = XDiffNorm;
		}
		
		//           std::cout << "ERROR Vector:\n"; error.Print();
		
		std::cout << "Convergence Order:\n";
		for(int j = 2; j < 11; j++)
		{
			std::cout << ( log(error(1,0)) - log(error(j,0)) )/( log(0.1)-log(j/10.) ) << "\n";
		}
		std::cout << "\nIf another kind of results are needed, edit the ConvTest Class on source code!\n";
		return;
	}
	
	else
	{
		std::cout << "Element don't fit in an option of Convergence Analys!\nSee ConvTest Class...\n";
		return;
	}
}
