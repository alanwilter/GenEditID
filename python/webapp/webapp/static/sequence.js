var infomap;

sortByName = function(a, b)
{
    return a.name.localeCompare(b.name);
}

labSelected = function()
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
    $("#newprojecttext").val('');
}

researcherSelected = function()
{
    var labselect = $('#lab_select');
    var researcherselect = $('#researcher_select');
    var projectselect = $('#project_select');

    projectselect.empty();
    $('#submitbutton').prop('disabled', true);
    setEnabled($('.projectsourceradio'), true);

    var labid = labselect.find(":selected").val();
    var researcherid = researcherselect.find(":selected").val();
    if (!!labid && !!researcherid)
    {
        var lab = infomap[labid];
        var researcher = lab.researchers[researcherid];

        var projectarray = $.map(researcher.projects, function(v) { return v; });
        
        if (projectarray.length > 0)
        {
            projectarray.sort(sortByName);
            
            for (var l = 0; l < projectarray.length; l++)
            {
                var project = projectarray[l];
                projectselect.append($("<option/>").val(project.id).text(project.name));
            }
            
            setEnabled($("#existingradio"), true);
            $("#existingradio").prop("checked", true).change();
        }
        else
        {
            $("#newradio").prop("checked", true).change();
            setEnabled($("#existingradio"), false);
        }
        $("#newprojecttext").val('');
        projectSourceChange();
    }
}

projectSourceChange = function()
{
    var useExisting = $('#existingradio').is(":checked")
    
    setEnabled($('#project_select'), useExisting);
    setEnabled($('#newprojecttext'), !useExisting);
}

projectSelected = function()
{
    var projectselect = $('#project_select');

    var projectId = projectselect.find(":selected").val();

    setEnabled($('#submitbutton'), !!projectId);
}

projectNameChange = function()
{
    var newProjectName = $("#newprojecttext").val();
    var hasName = newProjectName.length > 0;
    
    setEnabled($('#submitbutton'), hasName);
}

submitForSequencing = function()
{
    return true;
}

setEnabled = function(element, enabled)
{
    var prop = enabled ? false : 'disabled';
    element.prop('disabled', prop);
}

sequenceProjectReady = function()
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

    labselect.change(labSelected);
    researcherselect.change(researcherSelected);
    projectselect.change(projectSelected);
    $('.projectsourceradio').change(projectSourceChange);
    $('#submit_sequencing_form').submit(submitForSequencing);
    $('#newprojecttext').bind("change paste keyup", projectNameChange);

    setEnabled($("#newprojecttext"), false);
    setEnabled($('.projectsourceradio'), false);

    var selectedplate = $('#plateidspan').text();
    if (!!selectedplate)
    {
        $('#plate_select').val(selectedplate).prop('selected', true);
    }

    var selectedlab = $('#labidspan').text();
    if (!!selectedlab)
    {
        labselect.val(selectedlab).prop('selected', true).change();

        var selectedresearcher = $('#researcheridspan').text();
        if (!!selectedresearcher)
        {
            researcherselect.val(selectedresearcher).prop('selected', true).change();

            var projectsource = $('#projectsourcespan').text();
            
            if (projectsource == 'new')
            {
                $("#newradio").prop("checked", true);
                $('#newprojecttext').val($('#newprojectspan').text());
                projectSourceChange();
                projectNameChange();
            }
            else
            {
                $("#existingradio").prop("checked", true);
                
                var selectedproject = $('#projectidspan').text();
                if (!!selectedproject)
                {
                    projectselect.val(selectedproject).prop('selected', true).change();
                }
                projectSourceChange();
            }
        }
    }
}
