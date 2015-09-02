package bio.controllers;

import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;

@Controller
public class BrowseController
{
    @RequestMapping(value = "/browse", method = RequestMethod.GET)
    public String browse()
    {
        return "browse";
    }
}
