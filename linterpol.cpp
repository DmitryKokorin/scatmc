#include <cstdio>

#include "linterpol.h"

LinearInterpolation::LinearInterpolation(const Float min, const Float max, const size_t size) :
    m_min(min),
    m_max(max),
    m_resolution((max - min)/(size - 1)),
    m_data()
{
    m_data.resize(size);
}

Float LinearInterpolation::operator()(const Float& x) const
{
    if (x < m_min || x > m_max)
        return 0.;


	//locate index
	Float mu = (x - m_min) / m_resolution;
	size_t  idx = (size_t)(mu);

	if (idx + 1 == m_data.size())
	    return m_data[idx];

	mu  = mu - idx;

	return m_data[idx]*(1. - mu) + m_data[idx+1] * mu;
}


bool LinearInterpolation::load(const std::string& name)
{
	FILE* file = fopen(name.c_str(), "r");

	if (!file)
		return false;

	int size;

	int res = 0;

	res = fscanf(file, "%d\t%le\t%le", &size, &m_min, &m_max);

	if (EOF == res) {

		fclose(file);
		return false;
	}

	m_resolution = (m_max - m_min)/(size - 1);
	m_data.resize(size);


	Float value;

	for (int i = 0; i < size; ++i) {

		res = fscanf(file, "%le\t%le", &value, &m_data[i]);

		if (EOF == res && i != size) {

			fclose(file);
			return false;
		}
	}


	fclose(file);
	return true;
}

bool LinearInterpolation::save(const std::string& name)
{
	FILE* file = fopen(name.c_str(), "w");

	if (!file)
		return false;


	fprintf(file, "%d\t%.17e\t%.17e\n", (int)m_data.size(), m_min, m_max);

	for (size_t i = 0; i < m_data.size(); ++i) {

		fprintf(file, "%.17e\t%.17e\n", i*m_resolution, m_data[i]);
	}

	fclose(file);

	return true;
}

