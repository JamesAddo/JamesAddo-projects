package bio.services;

import bio.daos.CityDAO;
import bio.domains.City;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import javax.transaction.Transactional;
import java.util.List;

@Service
public class CityService
{
    @Autowired
    private CityDAO cityDAO;

    @Transactional
    public void insertCities(List<City> cities)
    {

    }
}
