using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PZ;

namespace BoundaryCross
{
    class Program
    {
        static void Main(string[] args)
        {

            Model.NormalDistribution();

            CalMatrix test = new CalMatrix();
            test.Values = new Double[,] { { 1,0,1 }, { 1, 1, 2 }, { 3, 4,2 } };
            test.CalDet();
            test.InverseMatrix();
            Double result = test.MatrixDet;
            Console.WriteLine("Please input the Direction and Velocity of the boundary: 4 5 6 0.01");
            String input = "4 5 6 0.01";//Console.ReadLine();
            String[] inputs = input.Split(' ');
            BoundaryInSpace test2 = new BoundaryInSpace(Double.Parse(inputs[0]), Double.Parse(inputs[1]), Double.Parse(inputs[2]), Double.Parse(inputs[3]),-30,-30,-30);
            test2.ResultCalwithTime();

            
        }
    }
    public class BoundaryInSpace
    {
        MovingFormation Tetrahedron = new MovingFormation();

        public void ResultCalwithTime()
        {
            //StreamWriter tw = File.CreateText(@"E:\MonteCarloForErrorAnalysis\Result0110\Result.txt");
            Double[] flags = new Double[] { 0, 0, 0, 0 };
            Double[] time = new Double[] { 0, 0, 0, 0 };
            //1,0,-1/sqrt(2)) , (-1,0,-1/sqrt(2)) , (0,1,1/sqrt(2)) , (0,-1,1/sqrt(2)
            Double[,] position = new Double[,] { { 1, 0, -1 / Math.Sqrt(2) },{ -1, 0, -1 / Math.Sqrt(2) }, { 0, 1, 1 / Math.Sqrt(2) }, { 0, -1, 1 / Math.Sqrt(2) } };
            double count = 4;
            for (t = 0; t < 10000; t+=0.001)
            {
                //tw.Write(t);
                
                for (int i = 0; i < 4; i++)
                {
                    Vector3D point = new Vector3D();
                    point = Tetrahedron.VertexPosition(t, i);
                    double tmpresult = this.CalAValueWithTimeAndPosition(point.x, point.y, point.z, t);
                    if (flags[i] == 0 && tmpresult == 1)
                    {
                        flags[i] = 1;
                        time[i] = t;
                       // position[i,0] = point.x;
                        //position[i, 1] = point.y;
                        //position[i, 2] = point.z;
                        count--;
                    }
                    //tw.Write(" ");
                    //tw.Write(tmpresult);
                }
                if (count <= 0)
                    break;
                //tw.WriteLine();
            }
            double t0 = (time[0] + time[1] + time[2] + time[3])/4;
            time[0] -= t0;
            time[1] -= t0;
            time[2] -= t0;
            time[3] -= t0;

            CalMatrix CalM = new CalMatrix();
            for(int i = 0;i<3;i++)
            {
                for(int j = 0;j < 3;j++)
                {
                    CalM.Values[i, j] = 0;
                    for(int k = 0;k<4;k++)
                    {
                        CalM.Values[i, j] += position[k, i] * position[k, j];
                    }
                }
            }
            Vector3D RightB = new Vector3D();
            RightB.x = 0;
            RightB.y = 0;
            RightB.z = 0;
            for (int i = 0;i<4;i++)
            {
                RightB.x += time[i] * position[i, 0];
                RightB.y += time[i] * position[i, 1];
                RightB.z += time[i] * position[i, 2];
            }
            CalM.CalDet();
            CalM.InverseMatrix();
            Vector3D mresult = Vector3D.VectorMultiMatrix(RightB, CalM);

            v = 1.0/ Math.Sqrt(mresult.x * mresult.x + mresult.y * mresult.y + mresult.z * mresult.z);
            mresult.x *= v;
            mresult.y *= v;
            mresult.z *= v;
            Console.WriteLine("The detected Vector of the Boundary is :");
            Console.Write(mresult.x);
            Console.Write("     ");
            Console.Write(mresult.y);
            Console.Write("     ");
            Console.Write(mresult.z);
            Console.WriteLine();
            Console.WriteLine("Detected velocity is : ");
            Console.Write(v);
            Vector3D yAxis = Vector3D.CrossProduct(mresult, new Vector3D());
            Vector3D zAxix = Vector3D.CrossProduct(yAxis, mresult);

            #region Single Parameter Chaning
            //for (int cout = 0; cout < 4; cout++)
            //{
            //    StreamWriter etw = File.CreateText(@"E:\MonteCarloForErrorAnalysis\Result0110\PositioningParameterChanging" + cout.ToString() + ".txt");
            //    for (int erroradd = 0; erroradd < 100; erroradd++)
            //    {
            //        Double[] TmpTime = new Double[4];
            //        for (int i = 0; i < 4; i++)
            //            TmpTime[i] = time[i];
            //        //TmpTime[cout] += 0.001 * erroradd;
            //        double[,] TmpPositions = new Double[,] { { 1, 0, -1 / Math.Sqrt(2) }, { -1, 0, -1 / Math.Sqrt(2) }, { 0, 1, 1 / Math.Sqrt(2) }, { 0, -1, 1 / Math.Sqrt(2) } };


            //        for (int j = 0; j < 3; j++)
            //            TmpPositions[cout, j] += 0.0001 * erroradd;

            //        CalMatrix TmpCalM = new CalMatrix();
            //        for (int i = 0; i < 3; i++)
            //        {
            //            for (int j = 0; j < 3; j++)
            //            {
            //                TmpCalM.Values[i, j] = 0;
            //                for (int k = 0; k < 4; k++)
            //                {
            //                    TmpCalM.Values[i, j] += TmpPositions[k, i] * TmpPositions[k, j];
            //                }
            //            }
            //        }
            //        Vector3D TmpRightB = new Vector3D();
            //        TmpRightB.x = 0;
            //        TmpRightB.y = 0;
            //        TmpRightB.z = 0;
            //        for (int i = 0; i < 4; i++)
            //        {
            //            TmpRightB.x += TmpTime[i] * TmpPositions[i, 0];
            //            TmpRightB.y += TmpTime[i] * TmpPositions[i, 1];
            //            TmpRightB.z += TmpTime[i] * TmpPositions[i, 2];
            //        }
            //        TmpCalM.CalDet();
            //        TmpCalM.InverseMatrix();
            //        Vector3D TmpResult = Vector3D.VectorMultiMatrix(TmpRightB, TmpCalM);
            //        double tmpv = 1.0 / Math.Sqrt(TmpResult.x * TmpResult.x + TmpResult.y * TmpResult.y + TmpResult.z * TmpResult.z);
            //        //TmpResult.x *= tmpv;
            //        //TmpResult.y *= tmpv;
            //        //TmpResult.z *= tmpv;
            //        double tmpAngle = Vector3D.AngleBetweenVector(TmpResult, mresult);
            //        double errorInVelocity = (tmpv - v) / v;
            //        etw.WriteLine((0.0001 * erroradd).ToString() + " " + TmpResult.x.ToString() + " " + TmpResult.y.ToString() + " " + TmpResult.z.ToString() + " " + errorInVelocity.ToString() + " " + tmpAngle.ToString());

            //    }
            //    etw.Flush();
            //    etw.Close();
            //}
            #endregion
            #region 2D map
            //for (int cout = 3; cout < 4; cout++)
            //{
            //    StreamWriter etw = File.CreateText(@"E:\MonteCarloForErrorAnalysis\Result0110\2DCover_V" /*+ cout.ToString() */+ ".txt");
            //    for (int xerroradd = 0; xerroradd < 200; xerroradd++)
            //    {
            //        for (int yerroradd = 0; yerroradd < 200; yerroradd++)
            //        {
            //            Double[] TmpTime = new Double[4];
            //            for (int i = 0; i < 4; i++)
            //                TmpTime[i] = time[i];
            //            TmpTime[0] = time[0] + 0.01 * xerroradd;
            //            double[,] TmpPositions = new Double[,] { { 1, 0, -1 / Math.Sqrt(2) }, { -1, 0, -1 / Math.Sqrt(2) }, { 0, 1, 1 / Math.Sqrt(2) }, { 0, -1, 1 / Math.Sqrt(2) } };


            //            for (int j = 0; j < 3; j++)
            //                TmpPositions[0, j] += 0.0001 * yerroradd;

            //            CalMatrix TmpCalM = new CalMatrix();
            //            for (int i = 0; i < 3; i++)
            //            {
            //                for (int j = 0; j < 3; j++)
            //                {
            //                    TmpCalM.Values[i, j] = 0;
            //                    for (int k = 0; k < 4; k++)
            //                    {
            //                        TmpCalM.Values[i, j] += TmpPositions[k, i] * TmpPositions[k, j];
            //                    }
            //                }
            //            }
            //            Vector3D TmpRightB = new Vector3D();
            //            TmpRightB.x = 0;
            //            TmpRightB.y = 0;
            //            TmpRightB.z = 0;
            //            for (int i = 0; i < 4; i++)
            //            {
            //                TmpRightB.x += TmpTime[i] * TmpPositions[i, 0];
            //                TmpRightB.y += TmpTime[i] * TmpPositions[i, 1];
            //                TmpRightB.z += TmpTime[i] * TmpPositions[i, 2];
            //            }
            //            TmpCalM.CalDet();
            //            TmpCalM.InverseMatrix();
            //            Vector3D TmpResult = Vector3D.VectorMultiMatrix(TmpRightB, TmpCalM);
            //            double tmpv = 1.0 / Math.Sqrt(TmpResult.x * TmpResult.x + TmpResult.y * TmpResult.y + TmpResult.z * TmpResult.z);
            //            //TmpResult.x *= tmpv;
            //            //TmpResult.y *= tmpv;
            //            //TmpResult.z *= tmpv;
            //            double tmpAngle = Vector3D.AngleBetweenVector(TmpResult, mresult);
            //            //etw.WriteLine((0.001 * xerroradd).ToString() + " " + (0.0001 * yerroradd).ToString() + " " /*+ TmpResult.x.ToString() + " " + TmpResult.y.ToString() + " " + TmpResult.z.ToString() + " "*/ + tmpAngle.ToString());
            //            etw.Write(/*(0.001 * xerroradd).ToString() + " " + (0.0001 * yerroradd).ToString() + " " *//*+ TmpResult.x.ToString() + " " + TmpResult.y.ToString() + " " + TmpResult.z.ToString() + " "*/ /*+*/ Math.Abs(v-tmpv).ToString() + " ");

            //        }
            //        etw.WriteLine();
            //    }
            //    etw.Flush();
            //    etw.Close();
            //}

            #endregion

// For Vector detection, the Limit-State-Function is assumed as G(x)=Angle|Detected - True| - 0.05
            #region MonteCarloCalculate
            //for(int MonteI = 0; MonteI<100; MonteI+=10)
            {

                StreamWriter tw = File.CreateText(@"E:\MonteCarloForErrorAnalysis\Result0110\MonteCarloAngleResult_All_Error_WRAMethodValidation" + /*MonteI.ToString()+*/ ".txt");
                //Random Error Generate
                GaussianRNG tmpRandom = new GaussianRNG();
                //.WriteLine("0 0 0 0 " + mresult.x.ToString() + " " + mresult.y.ToString() + " " + mresult.z.ToString() + " " + v.ToString());
                tw.WriteLine(mresult.x.ToString() + " " + mresult.y.ToString() + " " + mresult.z.ToString());
                Random rforV = new Random(unchecked((int)DateTime.Now.Ticks));
                for (int RepeatTimes = 0; RepeatTimes < 100000; RepeatTimes++)
                {

                    Double[] RandomErrorForTime1 = new Double[2];//Model.NormalDistribution();
                    Double[] RandomErrorForTime2 = new Double[2];//Model.NormalDistribution();
                    RandomErrorForTime1[0] = 0.1* 0.1 * tmpRandom.Next();
                    RandomErrorForTime2[0] = 0.1 * 0.3 * tmpRandom.Next();
                    RandomErrorForTime1[1] = 0.1 * 0.2 * tmpRandom.Next();
                    RandomErrorForTime2[1] = 0.1 * 0.4 * tmpRandom.Next();
                    Double[] TmpTime = new Double[4];
                    for (int i = 0; i < 4; i++)
                        TmpTime[i] = time[i];
                    //Add Timing error for Sat0 -Sat3

                    TmpTime[0] += /*MonteI **/ RandomErrorForTime1[0];
                    TmpTime[1] += /*MonteI **/ RandomErrorForTime1[1];
                    TmpTime[2] += /*MonteI **/ RandomErrorForTime2[0];
                    TmpTime[3] += /*MonteI **/ RandomErrorForTime2[1];

                    Vector3D TmpRightB = new Vector3D();
                    TmpRightB.x = 0;
                    TmpRightB.y = 0;
                    TmpRightB.z = 0;
                    #region position error
                    
                    double[,] TmpPositions = new Double[,] { { 1, 0, -1 / Math.Sqrt(2) }, { -1, 0, -1 / Math.Sqrt(2) }, { 0, 1, 1 / Math.Sqrt(2) }, { 0, -1, 1 / Math.Sqrt(2) } };
                    //Add Positioning error for Sat0 - Sat3
                    Double[] RandomErrorForP = new Double[4];
                    RandomErrorForP[0] = 0.004 * tmpRandom.Next();
                    RandomErrorForP[1] = 0.003 * tmpRandom.Next();
                    RandomErrorForP[2] = 0.002 * tmpRandom.Next();
                    RandomErrorForP[3] = 0.001 * tmpRandom.Next();
                    for (int i = 0; i < 4; i++)
                    {
                        #region Generate Random Vector
                        double u = rforV.Next(0,1);
                        double v = rforV.Next(0, 1);
                        double theta = 2 * Math.PI * u;
                        double phi = Math.Acos(2 * v - 1);
                        Vector3D randomv = new Vector3D();
                        randomv.x = Math.Sin(theta) * Math.Sin(phi);
                        randomv.y = Math.Cos(theta) * Math.Sin(phi);
                        randomv.z = Math.Cos(phi);
                        #endregion
                        TmpPositions[i, 0] += RandomErrorForP[i] * randomv.x;
                        TmpPositions[i, 1] += RandomErrorForP[i] * randomv.y;
                        TmpPositions[i, 2] += RandomErrorForP[i] * randomv.z;
                        //for (int j = 0; j < 3; j++)
                        //    TmpPositions[i, j] += 0.001/**MonteI*/ * tmpRandom.Next();
                    }
                    CalMatrix TmpCalM = new CalMatrix();
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            TmpCalM.Values[i, j] = 0;
                            for (int k = 0; k < 4; k++)
                            {
                                TmpCalM.Values[i, j] += TmpPositions[k, i] * TmpPositions[k, j];
                            }
                        }
                    }
                    TmpCalM.CalDet();
                    TmpCalM.InverseMatrix();

                    #endregion
                    for (int i = 0; i < 4; i++)
                    {
                        TmpRightB.x += TmpTime[i] * TmpPositions[i, 0];
                        TmpRightB.y += TmpTime[i] * TmpPositions[i, 1];
                        TmpRightB.z += TmpTime[i] * TmpPositions[i, 2];
                    }
                    //CalM.CalDet();
                    //CalM.InverseMatrix();
                    Vector3D TmpResult = Vector3D.VectorMultiMatrix(TmpRightB, TmpCalM);

                    double tmpv = 1.0 / Math.Sqrt(TmpResult.x * TmpResult.x + TmpResult.y * TmpResult.y + TmpResult.z * TmpResult.z);
                    //TmpResult.x *= tmpv;
                    //TmpResult.y *= tmpv;
                    //TmpResult.z *= tmpv;
                    double tmpAngle = Vector3D.AngleBetweenVector(TmpResult, mresult);

                    double Aresult = tmpAngle - 0.05;
                    //tw.WriteLine(RandomErrorForTime1[0].ToString() + " " + RandomErrorForTime1[1].ToString() + " " + RandomErrorForTime2[0].ToString() + " " + RandomErrorForTime2[1].ToString() + " " + TmpResult.x.ToString() + " " + TmpResult.y.ToString() + " " + TmpResult.z.ToString() + " " + (tmpv - v).ToString() + " " + tmpAngle.ToString());
                    //tw.WriteLine(TmpResult.x.ToString() + " " + TmpResult.y.ToString() + " " + TmpResult.z.ToString() + " " + (tmpv - v).ToString() + " " + tmpAngle.ToString());

                    tw.WriteLine(RandomErrorForTime1[0].ToString() + " " + RandomErrorForTime1[1].ToString() + " " + RandomErrorForTime2[0].ToString() + " " + RandomErrorForTime2[1].ToString() + " " 
                        + RandomErrorForP[0].ToString() +  " " + RandomErrorForP[1].ToString() + " " + RandomErrorForP[2].ToString() + " " + RandomErrorForP[3].ToString() + " "
                        + Aresult.ToString());
                    //+ (tmpv - v).ToString());

                    //Double[] RandomErrorForPos1 = Model.NormalDistribution();
                    //Double[] RandomErrorForPos2 = Model.NormalDistribution();

                }
                tw.Close();
            }
            #endregion
            double canubelieveit = 0;
            //tw.Close();
        }
        double SpaceXmin;// = 0;
        double SpaceXmax;// = 10000;
        double SpaceYmin;// = 0;
        double SpaceYmax;// = 10000;
        double SpaceZmin;// = 0;
        double SpaceZmax;// = 10000;

        //Assume that there is a plane boundary moving at volecity v in the direction of n; Thickness is 5;
        double nx;// = 1;
        double ny;// = 1;
        double nz;// = 1;
        double nThickness;// = 5;
        double v;// = 0.1;
        double nLength;
        //At time t0, the boundary is at the origin point 0,0,0;
        double nxt0; //= 0;
        double nyt0; //= 0;
        double nzt0; // 0;
        double t;
        public BoundaryInSpace()
        {
            SpaceXmin = 0;
            SpaceXmax = 10000;
            SpaceYmin = 0;
            SpaceYmax = 10000;
            SpaceZmin = 0;
            SpaceZmax = 10000;

            nx = 1/Math.Sqrt(3);
            ny = 0;// 1 / Math.Sqrt(3);
            nz = 0;// 1 / Math.Sqrt(3);

            
            v = 0.01;

            nxt0 = -30;
            nyt0 = -30;
            nzt0 = -30;
            Init();
        }
        public BoundaryInSpace(double inx, double iny, double inz, double iv, double inxt0, double inyt0, double inzt0)
        {
            nx = inx / Math.Sqrt(inx*inx + iny*iny + inz*inz);
            ny = iny / Math.Sqrt(inx * inx + iny * iny + inz * inz);
            nz = inz / Math.Sqrt(inx * inx + iny * iny + inz * inz);
            v = iv;
            nxt0 = inxt0;
            nyt0 = inyt0;
            nzt0 = inzt0;
            Init();

        }
        void Init()
        {
            nThickness = 5;
            nLength = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            this.nx = nx / nLength;
            this.ny = ny / nLength;
            this.nz = nz / nLength;
        }
        public double CalAValueWithTimeAndPosition(double inputx, double inputy, double inputz, double t)
        {
            double sq = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            double D1 = Math.Abs(nx * inputx + ny * inputy + nz * inputz - nxt0 * nx - nyt0 * ny - nzt0 * nz - v*t) / sq;
            double D2 = Math.Abs(nx * inputx + ny * inputy + nz * inputz - nxt0 * nx - nyt0 * ny - nzt0 * nz - 5 - v * t) / sq;
            if (D1 <= nThickness && D2 <= nThickness)
                return 1;
            else
                return 0;
        }
    }

    public class Vector3D
    {
        public double x;
        public double y;
        public double z;
        public Vector3D()
        {
            x = 1;
            y = 0;
            z = 0;
        }
        public Double Length()
        {
            return Math.Sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
        }
        public static Vector3D VectorMultiMatrix(Vector3D B, CalMatrix IA)
        {
            Vector3D result = new Vector3D();
            result.x = IA.CreatedM[0, 0] * B.x + IA.CreatedM[0, 1] * B.y + IA.CreatedM[0, 2] * B.z;
            result.y = IA.CreatedM[1, 0] * B.x + IA.CreatedM[1, 1] * B.y + IA.CreatedM[1, 2] * B.z;
            result.z = IA.CreatedM[2, 0] * B.x + IA.CreatedM[2, 1] * B.y + IA.CreatedM[2, 2] * B.z;
            return result;
        }
        public static Double AngleBetweenVector(Vector3D A, Vector3D B)
        {
            return (180.0 / Math.PI)*Math.Acos((A.x*B.x+A.y*B.y+A.z*B.z)/(A.Length()*B.Length()));
        }
        public static Vector3D CrossProduct(Vector3D A, Vector3D B)
        {
            Vector3D result = new Vector3D();
            result.x = A.y * B.z - A.z * B.y;
            result.y = A.z * B.x - A.x * B.z;
            result.z = A.x * B.y - A.y * B.x;
            return result;
        }
    }
    public class MovingFormation
    {
        Double CenterX;
        Double CenterY;
        Double CenterZ;

        Double Volecity;
        Double Mnx;
        Double Mny;
        Double Mnz;

        //The vector of four vertexs are (1,0,-1/sqrt(2)) , (-1,0,-1/sqrt(2)) , (0,1,1/sqrt(2)) , (0,-1,1/sqrt(2))
        Double FormationScale;
        //Double P1x = 1;
        //Double P1y = 0;
        //Double P1z = -1 / Math.Sqrt(2);
        //Double P2x = -1;
        //Double P2y = 0;
        //Double P2z = -1 / Math.Sqrt(2);
        //Double P3x = 0;
        //Double P3y = 1;
        //Double P3z = 1 / Math.Sqrt(2);
        //Double P4x = 0;
        //Double P4y = -1;
        //Double P4z = 1 / Math.Sqrt(2);

        List<Double> Px = new List<double>();
        List<Double> Py = new List<double>();
        List<Double> Pz = new List<double>();
        public MovingFormation()
        {
            FormationScale = 1;
            Volecity = 0.1;
            Mnx = 0;
            Mny = 0;
            Mnz = 0;

            CenterX = 0;
            CenterY = 0;
            CenterZ = 0;
            Init();
        }
        void Init()
        {
            Px.Add(1);
            Px.Add(-1);
            Px.Add(0);
            Px.Add(0);

            Py.Add(0);
            Py.Add(0);
            Py.Add(1);
            Py.Add(-1);

            Pz.Add(-1 / Math.Sqrt(2));
            Pz.Add(-1 / Math.Sqrt(2));
            Pz.Add(1 / Math.Sqrt(2));
            Pz.Add(1 / Math.Sqrt(2));
        }

        public Vector3D VertexPosition(double t, int VertexIndex)
        {
            CenterX = Mnx * Volecity * t;
            CenterY = Mny * Volecity * t;
            CenterZ = Mnz * Volecity * t;
            Vector3D result = new Vector3D();
            result.x = CenterX + Px[VertexIndex];
            result.y = CenterY + Py[VertexIndex];
            result.z = CenterZ + Pz[VertexIndex];
            return result;
        }

    }
    public class Matrix3D
    {
        Vector3D XAxis = new Vector3D();
        Vector3D YAxis = new Vector3D();
        Vector3D ZAxis = new Vector3D();
        public Matrix3D()
        {
            YAxis.x = 0;
            YAxis.y = 1;
            YAxis.z = 0;
            ZAxis.x = 0;
            ZAxis.y = 0;
            ZAxis.z = 1;
        }
        public Matrix3D(Vector3D X,Vector3D Y,Vector3D Z)
        {
            this.XAxis = X;
            this.YAxis = Y;
            this.ZAxis = Z;
        }
        public Vector3D VectorFromOrigin(Vector3D input)
        {
            Vector3D result = new Vector3D();
            result.x = this.XAxis.x * input.x + this.YAxis.x * input.x + this.ZAxis.x * input.x;
            result.y = this.XAxis.y * input.x + this.YAxis.y * input.x + this.ZAxis.x * input.x;

            result.z = this.XAxis.x * input.x + this.YAxis.x * input.x + this.ZAxis.x * input.x;
            return result;
        }
    }
    public class CalMatrix
    {
        int n = 3;
        public Double[,] Values;
        public Double MatrixDet;
        public Double[,] CreatedM;
        Double[,] MidM;
        public CalMatrix()
        {
            n = 3;
            Init();
        }
        public CalMatrix(int inputn)
        {
            n = inputn;
            Init();
        }
        void Init()
        {
            Values = new Double[n, n];
            CreatedM = new Double[n, n];
            MidM = new Double[n - 1, n - 1];
        }
        public void CalDet()
        {
            int r, c, m;
            int lop = 0;
            double result = 0;
            double mid = 1;
            if (n != 1)
            {
                lop = (n == 2) ? 1 : n;            //控制求和循环次数,若为2阶，则循环1次，否则为n次  
                for (m = 0; m < lop; m++)
                {
                    mid = 1;            //顺序求和, 主对角线元素相乘之和  
                    for (r = 0, c = m; r < n; r++, c++)
                    {
                        mid = mid * Values[r, c % n];//(*(p + r * n + c % n));
                    }
                    result += mid;
                }
                for (m = 0; m < lop; m++)
                {
                    mid = 1;            //逆序相减, 减去次对角线元素乘积  
                    for (r = 0, c = n - 1 - m + n; r < n; r++, c--)
                    {
                        mid = mid * Values[r, c % n];//(*(p + r * n + c % n));
                    }
                    result -= mid;
                }
            }
            else
                result = Values[0,0];
            MatrixDet =  result;
        }

        public Double CalMidDet()
        {
            int midn = n - 1;
            int r, c, m;
            int lop = 0;
            double result = 0;
            double mid = 1;
            if (midn != 1)
            {
                lop = (midn == 2) ? 1 : midn;            //控制求和循环次数,若为2阶，则循环1次，否则为n次  
                for (m = 0; m < lop; m++)
                {
                    mid = 1;            //顺序求和, 主对角线元素相乘之和  
                    for (r = 0, c = m; r < midn; r++, c++)
                    {
                        mid = mid * MidM[r, c % midn];//(*(p + r * n + c % n));
                    }
                    result += mid;
                }
                for (m = 0; m < lop; m++)
                {
                    mid = 1;            //逆序相减, 减去次对角线元素乘积  
                    for (r = 0, c = midn - 1 - m + midn; r < midn; r++, c--)
                    {
                        mid = mid * MidM[r, c % midn];//(*(p + r * n + c % n));
                    }
                    result -= mid;
                }
            }
            else
                result = Values[0, 0];
            return result;
        }

        public void InverseMatrix()
        {
            for (int i = 0; i < n; i++)    //求逆矩阵  
            {
                for (int j = 0; j < n; j++)
                {
                    CreatedM[j,i] = Creat_M(i, j,n) / MatrixDet;
                }
            }
        }
        Double Creat_M( int m, int n,int k)
        {
            int len;
            int i, j;
            double mid_result = 0;
            int sign = 1;
            //float* p_creat, *p_mid;
            len = (k - 1) * (k - 1);            //k阶矩阵的代数余之式为k-1阶矩阵  
            //p_creat = (float*)calloc(len, sizeof(float)); //分配内存单元  
            //p_mid = p_creat;
            for (i = 0; i < k; i++)
            {
                for (j = 0; j < k; j++)
                {
                    if (i != m && j != n) //将除第i行和第j列外的所有元素存储到以p_mid为首地址的内存单元  
                    {
                        MidM[i > m ? i - 1 : i, j > n ? j - 1 : j] = Values[i,j];
                        //*p_mid++ = *(p + i * k + j);
                    }
                }
            }
            sign = (m + n) % 2 == 0 ? 1 : -1;    //代数余之式前面的正、负号  
            mid_result = (Double)sign * CalMidDet();
            //free(p_creat);
            return mid_result;
        }
    }
}
