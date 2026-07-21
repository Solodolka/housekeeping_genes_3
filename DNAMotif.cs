using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
//using System.Threading.Tasks;

namespace DNAMotif_01
{
    public struct Motif
    {
        public int Position;
        public double Weight;

        public Motif(int pos, double w)
        {
            Position = pos;
            Weight = w;
        }
    }


    public class DNAMotif
    {
        double[] matrix;
        int rowLength;
        int rows = 4;
        double threshold = 7.0;

        private Dictionary<char, int> nucToIdx = 
            new Dictionary<char, int> { { 'A', 0 }, { 'C', 1 }, { 'G', 2 }, { 'T', 3 } };

        public void SetMatrix(double[] m)
        {
            matrix = m;
            rowLength = matrix.Length / rows;

            Console.WriteLine($"Matrix {rows} x {rowLength}:");

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < rowLength; j++)
                {
                    // Формула индекса для плоского массива: [номер_строки * ширина + номер_столбца]
                    double value = matrix[i * rowLength + j];
                    Console.Write($"{value:F2}\t"); // Вывод с 2 знаками после запятой и табуляцией
                }
                Console.WriteLine(); // Переход на новую строку после завершения текущей
            }
        }

        public void SetThreshold(double t)
        {
            threshold = t;
        }

        public List<Motif> FindMotifs(string record)
        {
            var results = new List<Motif>();

            // Перебираем последовательность скользящим окном
            for (int i = 0; i < record.Length - rowLength + 1; i++)
            {
                double score = 0;
                // Вычисляем вес (score) для текущего окна
                for (int j = 0; j < rowLength; j++)
                {
                    char nuc = record[i + j];
                    if (nucToIdx.ContainsKey(nuc))
                    {
                        int idx = nucToIdx[nuc];
                        score += matrix[idx * rowLength + j];
                    }
                }

                score /= rowLength;

                if (score >= threshold)
                {
                    results.Add(new Motif(i, score));
                }
            }
            return results;

        }

        public List<Motif> FindMotifsTATA1999(string record)
        {
            var results = new List<Motif>();

            // Перебираем последовательность скользящим окном
            for (int i = 0; i < record.Length - rowLength + 1; i++)
            {
                double score = 0;
                // Вычисляем вес (score) для текущего окна
                for (int j = 0; j < rowLength; j++)
                {
                    char nuc = record[i + j];
                    if (j == 1 && nuc == 'G') score -= 100;
                    if (j == 1 && nuc == 'C') score -= 100;
                    if (j == 2 && nuc == 'A') score -= 100;
                    if (j == 3 && nuc == 'G') score -= 100;
                    if (j == 4 && nuc == 'C') score -= 100;
                    if (nucToIdx.ContainsKey(nuc))
                    {
                        int idx = nucToIdx[nuc];
                        score += matrix[idx * rowLength + j];
                    }
                }

                score /= rowLength;

                if (score >= threshold)
                {
                    results.Add(new Motif(i, score));
                }
            }
            return results;

        }
    }
}
