package bio.daos;

import org.hibernate.SessionFactory;
import bio.domains.BioActiveConstituent;
import org.springframework.stereotype.Repository;

import javax.annotation.Resource;

@Repository
public class BioActiveConstituentDAOHibernateImpl implements BioActiveConstituentDAO
{
    @Resource
    private SessionFactory sessionFactory;

    @Override
    public void persistBioActiveConstituent(BioActiveConstituent activeConstituent)
    {
        sessionFactory.getCurrentSession().save(activeConstituent);
    }
}
