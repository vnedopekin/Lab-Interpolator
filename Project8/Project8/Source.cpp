#include <iostream>
#include <fstream>
#include <cmath>
#include<array>

template<typename xType, typename yType, unsigned int N>

class NewtonInterpolator 
{
    std::array<xType, N> X;
    std::array<yType, N> Y;
public:
    NewtonInterpolator(const std::array<xType, N>& points, const std::array<yType, N>& values) noexcept {

        X = points;
        Y = values;


        // вычисление разделённых разностей
        for (int i = 0; i < N - 1; i++) 
        {
            for (int j = N - 1; j > i; j--) 
            {
                Y[j] = (Y[j] - Y[j - 1]) / (X[j] - X[j - 1 - i]);
            }

        }

    }

    //Вычисление полинома Ньютона по схеме Горнера
    yType interpolate(const xType& x) const noexcept {
        yType sum = 0;
        for (int i = N - 1; i > 0; i--) 
        {
            sum += Y[i];
            sum *= (x - X[i - 1]);
        }
        sum += Y[0];
        return sum;
    }

};

int main() {

    const unsigned int N = 4; // количество узлов
    std::ofstream fout;
    fout.open("4points.txt"); 


    for (int k = 1; k > (-5); k--) {

        std::array<double, N> X;
        std::array<double, N> Y;

        // заполнение массива координат x, y
        for (int i = 0; i < N; i++) {
            X[i] = (pow(2, k) * i) / (N - 1);
            Y[i] = exp(X[i]);
        }

        //интерполянт на равномерном распределении узлов
        NewtonInterpolator<double, double, N> M(X, Y);


        //максимальная и текущая ошибки для равномерного распределения узлов
        double error_even = 0;
        double current = 0;


        //ошибка как максимальная разница в 1000 точках, равномерно распределённых для интерполянта с равномерным распределением
        for (int i = 0; i < 1000; i++) {
            current = abs(exp(i / 999. * pow(2, k)) - M.interpolate(i / 999. * pow(2, k)));
            if (error_even < current) error_even = current;
        }



        // Для нулей полинома Чебышева
        std::array<double, N> Xch;
        std::array<double, N> Ych;


        //заполнение массивов для координат x, y
        for (int i = 0; i < N; i++) {
            Xch[i] = pow(2, k - 1) + pow(2, k - 1) * cos((2 * i + 1) * 3.14159265358979323846 / (2 * N));
            Ych[i] = exp(Xch[i]);
        }

        NewtonInterpolator<double, double, N> Mch(Xch, Ych);

        //максимальная и текущая ошибка для Чебышевского распределения
        double error_cheb = 0;
        double current_cheb = 0;



        //ошибка как максимальная разница в 1000 точках, равномерно распределенных для Чебышевского интерполянта
        for (int i = 0; i < 1000; i++) {
            current_cheb = abs(exp(i / 999. * pow(2, k)) - Mch.interpolate(i / 999. * pow(2, k)));
            if (error_cheb < current_cheb) error_cheb = current_cheb;
        }

        fout << error_even << " " << error_cheb << " " << std::endl;


    }
    fout.close();

    return 0;
}