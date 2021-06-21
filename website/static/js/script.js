// Code goes here

angular.module('bootstrapTabsDemo', ['ui.route']).
config(function ($routeProvider) {
  $routeProvider.when('/tab1', {
    template: '<div>Tab1</div>'
  }).when('/tab2', {
    template: '<div>Tab2</div>'
  });
});
