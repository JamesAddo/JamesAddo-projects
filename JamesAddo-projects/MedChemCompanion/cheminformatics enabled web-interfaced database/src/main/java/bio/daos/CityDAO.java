package bio.daos;

import bio.domains.City;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.stereotype.Repository;

import javax.annotation.Resource;

@Repository
public class CityDAO
{
    @Resource
    private JdbcTemplate jdbcTemplate;

    public void insertCity(City city)
    {
        jdbcTemplate.update(
            "INSERT INTO city VALUE(?, ? , ?)",
            new Object[] {
                city.getName(),
                city.getZipcode()
            }
        );
    }
}
