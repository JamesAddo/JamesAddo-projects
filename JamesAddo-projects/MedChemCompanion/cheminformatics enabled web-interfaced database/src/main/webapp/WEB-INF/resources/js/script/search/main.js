require.config({
    baseUrl: 'resources/js',
    paths: {
        'jquery': 'lib/jquery/jquery-1.11.0',
        'angular': 'lib/angular/angular',
        'angular-route': 'lib/angular/angular-route',
        'chemaxon-util': 'lib/chemaxon/util',
        'chemaxon-webservices': 'lib/chemaxon/webservices'
    },
    shim: {
        'angular': {
            deps: ['jquery'],
            exports: 'angular'
        },
        'angular-route': {
            deps: ['angular'],
            exports: 'angular-route'
        },
        'chemaxon-util': {
            deps: ['jquery'],
            export: 'chemaxon-util'
        },
        'chemaxon-webservices': {
            export: 'chemaxon-util'
        }
    }
});

require(['jquery', 'angular', 'angular-route', 'chemaxon-util', 'chemaxon-webservices'], function($, A) {
    'use strict';

    $(document).ready(function() {
        A.module('app', ['ngRoute']);

        getMarvinPromise("#sketch").done(function (sketcherInstance) {
            $("#exportBtn").on("click", function handleGetSmilesButton () {
                var getMolConvertURL = function () {
                    var ws = getDefaultServices();
                    return ws['molconvertws'];
                };

                var	data = JSON.stringify({
                    "structure" : sketcherInstance.exportAsMrv(),
                    "inputFormat" : "mrv",
                    "parameters" : "smiles"
                });

                $.ajax({
                    "url": getMolConvertURL()
                    ,"type": "POST"
                    ,"dataType": "json"
                    ,"contentType": "application/json"
                    ,"data": data
                }).done(function (data, textStatus, jqXHR) {
                    $("#searchBar").val(data['structure']);
                }).fail(function() {
                    console.log("Something wrong!");
                });
            });
        }).fail(function () {
            alert("Cannot retrieve sketcher instance from iframe");
        });

        $("#searchBtn").click(function() {
            $.ajax({
                type: "POST",
                url: "/search",
                data: { "query": $("#searchBar").val() },
                success: function(data) {
                    console.log(data);
                }
            });
        });
    });
});
