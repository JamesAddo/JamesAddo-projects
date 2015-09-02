package bio.domains;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.Table;
import java.io.Serializable;

@Entity
@Table(name="BIOACTIVE_CONSTITUENTS")
public class BioActiveConstituent implements Serializable
{
    @Id
    @Column(name="CD_ID")
    private Long id;

    @Column(name="ENTRY_NAME")
    private String name;

    @Column(name="CAS_REGISTRY_NUMBER")
    private String registryNumber;

    @Column(name="CD_SMILES")
    private String smiles;

    public Long getId()
    {
        return id;
    }

    public void setId(Long id)
    {
        this.id = id;
    }

    public String getName()
    {
        return name;
    }

    public void setName(String name)
    {
        this.name = name;
    }

    public String getRegistryNumber()
    {
        return registryNumber;
    }

    public void setRegistryNumber(String registryNumber)
    {
        this.registryNumber = registryNumber;
    }

    public void setSmiles(String smiles)
    {
        this.smiles = smiles;
    }

    public String getSmiles()
    {
        return smiles;
    }
}
