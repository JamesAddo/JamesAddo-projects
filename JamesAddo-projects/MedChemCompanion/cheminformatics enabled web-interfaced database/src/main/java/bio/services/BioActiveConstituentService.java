package bio.services;

import bio.daos.BioActiveConstituentDAO;
import bio.domains.BioActiveConstituent;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;

@Service
public class BioActiveConstituentService
{
    @Autowired
    private BioActiveConstituentDAO activeConstituentDAO;

    @Transactional
    public void saveActiveConstituent(BioActiveConstituent activeConstituent)
    {
        activeConstituentDAO.persistBioActiveConstituent(activeConstituent);
    }
}
