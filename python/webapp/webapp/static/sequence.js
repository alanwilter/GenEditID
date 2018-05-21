var infomap;

sortByName = function(a, b)
{
    return a.name.localeCompare(b.name);
}

labselected = function()
{
    var labselect = $('#lab_select');
    var researcherselect = $('#researcher_select');
    var projectselect = $('#project_select');

    researcherselect.empty();
    projectselect.empty();
    $('#submitbutton').prop('disabled', true);

    var labid = labselect.find(":selected").val();
    if (!!labid)
    {
        var lab = infomap[labid];

        var researcherarray = $.map(lab.researchers, function(v) { return v; });
        researcherarray.sort(sortByName);
        
        for (var l = 0; l < researcherarray.length; l++)
        {
            var researcher = researcherarray[l];
            researcherselect.append($("<option/>").val(researcher.id).text(researcher.name));
        }
    }
    
    researcherselect.prop('disabled', false);
    projectselect.prop('disabled', 'disabled');
}

researcherselected = function()
{
    var labselect = $('#lab_select');
    var researcherselect = $('#researcher_select');
    var projectselect = $('#project_select');

    projectselect.empty();
    $('#submitbutton').prop('disabled', true);

    var labid = labselect.find(":selected").val();
    var researcherid = researcherselect.find(":selected").val();
    if (!!labid && !!researcherid)
    {
        var lab = infomap[labid];
        var researcher = lab.researchers[researcherid];

        var projectarray = $.map(researcher.projects, function(v) { return v; });
        projectarray.sort(sortByName);
        
        for (var l = 0; l < projectarray.length; l++)
        {
            var project = projectarray[l];
            projectselect.append($("<option/>").val(project.id).text(project.name));
        }
    }
    
    projectselect.prop('disabled', false);
}

projectselected = function()
{
    var projectselect = $('#project_select');

    var projectid = projectselect.find(":selected").val();

    var prop = !!projectid ? false : 'disabled';
    $('#submitbutton').prop('disabled', prop);
}

submitting = function()
{
    alert("Submitting.");
    return true;
}

sequenceprojectready = function()
{
    infomap = JSON.parse($('#jsonmap').text());

    var labselect = $('#lab_select');
    var researcherselect = $('#researcher_select');
    var projectselect = $('#project_select');

    var labarray = $.map(infomap, function(v) { return v; });
    labarray.sort(sortByName);
    
    for (var l = 0; l < labarray.length; l++)
    {
        var lab = labarray[l];
        labselect.append($("<option/>").val(lab.id).text(lab.name));
    }

    labselect.change(labselected);
    researcherselect.change(researcherselected);
    projectselect.change(projectselected);
    $('#submit_sequencing_form').submit(submitting);

    var selectedlab = $('#labidspan').text();
    if (!!selectedlab)
    {
        labselect.val(selectedlab).prop('selected', true).change();

        var selectedresearcher = $('#researcheridspan').text();
        if (!!selectedresearcher)
        {
            researcherselect.val(selectedresearcher).prop('selected', true).change();

            var selectedproject = $('#projectidspan').text();
            if (!!selectedresearcher)
            {
                projectselect.val(selectedproject).prop('selected', true).change();
            }
        }
    }
}
