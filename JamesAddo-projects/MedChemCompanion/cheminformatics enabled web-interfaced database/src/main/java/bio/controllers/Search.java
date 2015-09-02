package bio.controllers;

import com.google.gson.Gson;
import bio.domains.BioActiveConstituent;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;

@Controller
@RequestMapping("/search")
public class Search
{
    @RequestMapping(method = RequestMethod.GET)
    public String get()
    {
        return "search";
    }

    @ResponseBody
    @RequestMapping(method = RequestMethod.POST)
    public String search(@RequestParam(value = "query") String query)
    {
        Gson g = new Gson();

        BioActiveConstituent result = new BioActiveConstituent();
        result.setId(1L);
        result.setName(query);
        result.setRegistryNumber("1");
        result.setSmiles("S1");

        return g.toJson(result);
    }
}
