using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
namespace PZ
{
    public class GaussianRNG
    {
        int iset;
        double gset;
        Random r1, r2;

        public GaussianRNG()
        {
            r1 = new Random(unchecked((int)DateTime.Now.Ticks));
            r2 = new Random(~unchecked((int)DateTime.Now.Ticks));
            iset = 0;
        }

        public double Next()
        {
            double fac, rsq, v1, v2;
            if (iset == 0)
            {
                do
                {
                    v1 = 2.0 * r1.NextDouble() - 1.0;
                    v2 = 2.0 * r2.NextDouble() - 1.0;
                    rsq = v1 * v1 + v2 * v2;
                } while (rsq >= 1.0 || rsq == 0.0);

                fac = Math.Sqrt(-2.0 * Math.Log(rsq) / rsq);
                gset = v1 * fac;
                iset = 1;
                return v2 * fac;
            }
            else
            {
                iset = 0;
                return gset;
            }
        }
    }
    class Model
    {
        /// <summary>
        /// 正态分布随机数
        /// </summary>
        const int N = 100;
        const int MAX = 50;
        const double MIN = 0.1;
        const int MIU = 40;
        const int SIGMA = 1;
        static Random aa = new Random((int)(DateTime.Now.Ticks / 10000));
        internal static double[] NormalDistribution()
        {
            Random rand = new Random((int)(DateTime.Now.Ticks / 100));
            double[] y;
            double u1, u2, v1 = 0, v2 = 0, s = 0, z1 = 0, z2 = 0;
            while (s > 1 || s == 0)
            {
                u1 = rand.NextDouble();
                u2 = rand.NextDouble();
                v1 = 2 * u1 - 1;
                v2 = 2 * u2 - 1;
                s = v1 * v1 + v2 * v2;
            }
            z1 = Math.Sqrt(-2 * Math.Log(s) / s) * v1;
            z2 = Math.Sqrt(-2 * Math.Log(s) / s) * v2;
            y = new double[] { z1, z2 };
            return y; //返回两个服从正态分布N(0,1)的随机数z0 和 z1
        }
        public double AverageRandom(double min, double max)//产生(min,max)之间均匀分布的随机数
        {
            int MINnteger = (int)(min * 10000);
            int MAXnteger = (int)(max * 10000);
            int resultInteger = aa.Next(MINnteger, MAXnteger);
            return resultInteger / 10000.0;
        }
        public double Normal(double x, double miu, double sigma) //正态分布概率密度函数
        {
            return 1.0 / (x * Math.Sqrt(2 * Math.PI) * sigma) * Math.Exp(-1 * (Math.Log(x) - miu) * (Math.Log(x) - miu) / (2 * sigma * sigma));
        }
        public double Random_Normal(double miu, double sigma, double min, double max)//产生正态分布随机数
        {
            double x;
            double dScope;
            double y;
            do
            {
                x = AverageRandom(min, max);
                y = Normal(x, miu, sigma);
                dScope = AverageRandom(0, Normal(miu, miu, sigma));
            } while (dScope > y);
            return x;
        }

        /// <summary>
        /// 指数分布随机数
        /// </summary>
        /// <param name="const_a"></param>
        /// <returns></returns>
        public double RandExp(double const_a)//const_a是指数分布的参数λ
        {
            Random rand = new Random(Guid.NewGuid().GetHashCode());
            double p;
            double temp;
            if (const_a != 0)
                temp = 1 / const_a;
            else
                throw new System.InvalidOperationException("除数不能为零！不能产生参数为零的指数分布！");
            double randres;
            while (true) //用于产生随机的密度，保证比参数λ小
            {
                p = rand.NextDouble();
                if (p < const_a)
                    break;
            }
            randres = -temp * Math.Log(temp * p, Math.E);
            return randres;
        }

        /// <summary>
        /// 负指数分布随机数产生
        /// </summary>
        /// <param name="lam">参数</param>
        /// <returns></returns>
        Random ran;
        public Model()
        {
            ran = new Random();
        }
        public double ngtIndex(double lam)
        {
            double dec = ran.NextDouble();
            while (dec == 0)
                dec = ran.NextDouble();
            return -Math.Log(dec) / lam;
        }

        /// <summary>
        /// 泊松分布产生
        /// </summary>
        /// <param name="lam">参数</param>
        /// <param name="time">时间</param>
        /// <returns></returns>
        public double poisson(double lam, double time)
        {
            int count = 0;
            while (true)
            {
                time -= ngtIndex(lam);
                if (time > 0)
                    count++;
                else
                    break;
            }
            return count;
        }

    }
}